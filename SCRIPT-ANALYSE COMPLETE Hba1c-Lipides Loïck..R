rm(list=ls(all=TRUE))
# ============================================================
# PLAN D'ANALYSE STATISTIQUE COMPLET — VERSION CORRIGÉE v2
# Profil glycémique HbA1c et corrélations lipidiques
# CUBT / CNHU-HKM · 1er juillet 2025 – 30 juin 2026 · n = 462
# ============================================================
# Corrections v2 :
#   [FIX-1] VIF polr : car::vif(m1a) direct (suppression workaround)
#   [FIX-2] Section B : vars_corr au lieu de vars_anal (hba1c exclus)
#   [FIX-3] Multicolinéarité CT+LDL : modèles a1/a2 séparés
#   [FIX-4] TG outliers : winsorisation + logistf alternatif
#   [FIX-5] ada_class coercion : factor() explicite
#   [FIX-6] aes_string deprecated : .data[[var]]
#   [FIX-7] nagelkerke_r2 : nobs(model) au lieu de model$data
# ============================================================


# ============================================================
# ÉTAPE 0 — PACKAGES ET IMPORT DES DONNÉES
# ============================================================

# install.packages(c(
#   "readxl","dplyr","tidyr","ggplot2","ggpubr",
#   "rstatix","nortest","DescTools","corrplot",
#   "MASS","brant","FSA","car","broom",
#   "pROC","ResourceSelection","nnet","logistf"
# ))

library(readxl); library(dplyr); library(tidyr)
library(ggplot2); library(ggpubr); library(rstatix)
library(nortest); library(DescTools); library(corrplot)
library(MASS); library(brant); library(FSA)
library(car); library(broom); library(pROC)
library(ResourceSelection); library(nnet); library(logistf)


# ============================================================
# ÉTAPE 1 — IMPORT ET RECODAGE
# ============================================================

df_raw <- read_excel(
  path      = "C:/Users/HP/Desktop/Ton dossier QGIS est ICI/Analyse/MON MEMOIRE/Base_CUBT_Robuste_2025-2026.xlsx",
  sheet     = "Donnees",

  col_names = TRUE
)
cat("Colonnes lues :", names(df_raw), "\n")
cat("Dimensions    :", dim(df_raw), "\n")

# --- 1.1 Renommage et typage -----------------------------------
df <- df_raw %>%
  dplyr::select(1:18) %>%
  setNames(c("id","doss","code","nom","age","sexe","date",
             "hba1c","ct","ldl","hdl","tg",
             "castelli","tg_hdl","ldl_hdl",
             "ada_class","cubt_class","notes")) %>%
  mutate(
    across(c(age,hba1c,ct,ldl,hdl,tg,castelli,tg_hdl,ldl_hdl),
           as.numeric),
    sexe = factor(sexe, levels=c("M","F"),
                  labels=c("Homme","Femme")),
    # [FIX-5] ada_class en factor dès l'import
    ada_class  = factor(ada_class,
                        levels=c("Normoglycémique","Pré-diabétique",
                                 "Diabétique","Hypoglycémie")),
    cubt_class = factor(cubt_class,
                        levels=c("Dans la plage","Hors plage")),
    # Variable ordonnée ADA
    ada_ord = factor(ada_class,
                     levels=c("Normoglycémique","Pré-diabétique",
                              "Diabétique"),
                     ordered=TRUE),
    # Variable binaire ADA
    ada_bin = factor(
      case_when(
        ada_class == "Normoglycémique"                     ~ "Non pathologique",
        ada_class %in% c("Pré-diabétique","Diabétique")   ~ "Pathologique",
        TRUE ~ NA_character_),
      levels=c("Non pathologique","Pathologique")),
    # Variable binaire CUBT
    cubt_bin = cubt_class,
    # Dyslipidémies dichotomisées
    ct_eleve       = as.integer(ct   > 2.00),
    ct_bas         = as.integer(ct   < 1.25),
    ldl_eleve      = as.integer(ldl  > 1.60),
    ldl_bas        = as.integer(ldl  < 1.00),
    hdl_eleve      = as.integer(hdl  > 0.70),
    hdl_bas        = as.integer(hdl  < 0.40),
    tg_eleve       = as.integer(tg   > 1.50),
    tg_bas         = as.integer(tg   < 0.50),
    castelli_patho = as.integer(
      (sexe=="Homme" & castelli>5.0) |
        (sexe=="Femme" & castelli>4.5)),
    tg_hdl_patho   = as.integer(tg_hdl  > 3.0),
    ldl_hdl_patho  = as.integer(ldl_hdl > 3.5)
  ) %>%
  filter(!is.na(hba1c))

# --- 1.2 Winsorisation TG au 99e percentile -------------------
# [FIX-4] : La valeur TG=44.00 g/l crée une séparation parfaite
# → on conserve les valeurs brutes dans df, on crée tg_wins pour la modélisation
p99_tg <- quantile(df$tg, 0.99, na.rm=TRUE)
cat(sprintf("\n[FIX-4] Seuil winsorisation TG (99e pct) : %.3f g/l\n", p99_tg))
cat(sprintf("Valeurs TG > %.3f : %d patients\n",
            p99_tg, sum(df$tg > p99_tg, na.rm=TRUE)))
df <- df %>%
  mutate(tg_wins = pmin(tg, p99_tg),
         tg_hdl_wins = round(tg_wins / hdl, 2))

# --- 1.3 Sous-ensembles analytiques ---------------------------
df_A     <- df %>%
  filter(ada_class %in% c("Normoglycémique","Pré-diabétique","Diabétique")) %>%
  droplevels()

df_A_bin <- df %>%
  filter(!is.na(ada_bin)) %>%
  droplevels()

df_B <- df    # Section B : n=462

cat("\n=== VÉRIFICATION DE LA BASE ===\n")
cat("n total           :", nrow(df),       "\n")
cat("n Section A (ord) :", nrow(df_A),     "\n")
cat("n Section A (bin) :", nrow(df_A_bin), "\n")
cat("n Section B       :", nrow(df_B),     "\n\n")
cat("Distribution ADA :\n");  print(table(df$ada_class))
cat("\nDistribution CUBT :\n"); print(table(df$cubt_class))
cat("\nValeurs manquantes :\n")
print(colSums(is.na(df[, c("hba1c","ct","ldl","hdl","tg",
                           "castelli","tg_hdl","ldl_hdl")])))


# ============================================================
# ÉTAPE 2 — STATISTIQUES DESCRIPTIVES
# ============================================================

describe_var <- function(x, var_name="") {
  tibble(
    Variable  = var_name,
    n         = sum(!is.na(x)),
    Moyenne   = round(mean(x,  na.rm=TRUE), 3),
    Ecart_type= round(sd(x,    na.rm=TRUE), 3),
    Mediane   = round(median(x,na.rm=TRUE), 3),
    Q1        = round(quantile(x, 0.25, na.rm=TRUE), 3),
    Q3        = round(quantile(x, 0.75, na.rm=TRUE), 3),
    Min       = round(min(x,   na.rm=TRUE), 3),
    Max       = round(max(x,   na.rm=TRUE), 3)
  )
}

vars_bio  <- c("age","hba1c","ct","ldl","hdl","tg",
               "castelli","tg_hdl","ldl_hdl")
# [FIX-2] vars_corr et vars_anal BIEN distincts
vars_corr <- c("ct","ldl","hdl","tg","castelli","tg_hdl","ldl_hdl")
labels_corr <- c("CT","LDL-c","HDL-c","TG","Castelli","TG/HDL","LDL/HDL")
vars_anal   <- c("hba1c",vars_corr)    # pour stats descriptives uniquement

cat("\n=== STATISTIQUES DESCRIPTIVES — POPULATION GLOBALE (n=462) ===\n")
desc_global <- bind_rows(lapply(vars_bio, function(v)
  describe_var(df[[v]], v)))
print(as.data.frame(desc_global), digits=3)

cat("\n=== STATISTIQUES PAR GROUPE ADA ===\n")
desc_ada <- df_A %>%
  group_by(ada_class) %>%
  summarise(n=n(),
            across(all_of(vars_corr),
                   list(moy=~round(mean(.x,na.rm=TRUE),3),
                        sd =~round(sd(.x,na.rm=TRUE),3),
                        med=~round(median(.x,na.rm=TRUE),3))),
            .groups="drop")
print(as.data.frame(desc_ada))

cat("\n=== PRÉVALENCES DES ANOMALIES (n=462) ===\n")
n_total <- nrow(df)
prev <- tibble(
  Parametre = c("CT élevé (>2.00)","CT bas (<1.25)",
                "LDL élevé (>1.60)","LDL bas (<1.00)",
                "HDL élevé (>0.70)","HDL bas (<0.40)",
                "TG élevé (>1.50)","TG bas (<0.50)",
                "Castelli pathologique","TG/HDL >3.0","LDL/HDL >3.5"),
  n = c(sum(df$ct_eleve),sum(df$ct_bas),
        sum(df$ldl_eleve),sum(df$ldl_bas),
        sum(df$hdl_eleve),sum(df$hdl_bas),
        sum(df$tg_eleve),sum(df$tg_bas),
        sum(df$castelli_patho),sum(df$tg_hdl_patho),sum(df$ldl_hdl_patho)),
  Pct = round(c(sum(df$ct_eleve),sum(df$ct_bas),
                sum(df$ldl_eleve),sum(df$ldl_bas),
                sum(df$hdl_eleve),sum(df$hdl_bas),
                sum(df$tg_eleve),sum(df$tg_bas),
                sum(df$castelli_patho),sum(df$tg_hdl_patho),
                sum(df$ldl_hdl_patho)) / n_total * 100, 1))
print(as.data.frame(prev))

cat("\n=== RÉPARTITION PAR TRANCHE D'ÂGE × SEXE ===\n")
df <- df %>%
  mutate(tranche_age = cut(age,
                           breaks=c(18,30,40,50),
                           labels=c("18-29 ans","30-39 ans","40-50 ans"),
                           include.lowest=TRUE))
print(table(df$tranche_age, df$sexe))

# --- 2.7 Graphiques -------------------------------------------
vars_plot   <- c("hba1c","ct","ldl","hdl","tg")
labels_plot <- c("HbA1c (%)","CT (g/l)","LDL-c (g/l)",
                 "HDL-c (g/l)","TG (g/l)")

# [FIX-6] aes_string() → .data[[var]]
# G1 : Histogrammes
p_hist_list <- lapply(seq_along(vars_plot), function(i) {
  v <- vars_plot[i]
  ggplot(df, aes(x=.data[[v]])) +
    geom_histogram(bins=25, fill="#12B3C4", color="#0D1F3A", alpha=0.8) +
    labs(title=labels_plot[i], x=labels_plot[i], y="Effectif") +
    theme_minimal(base_size=11) +
    theme(plot.title=element_text(hjust=0.5, face="bold"))
})
p_hist <- ggarrange(plotlist=p_hist_list, ncol=3, nrow=2)
ggsave("fig1_histogrammes.png", p_hist, width=14, height=9, dpi=300)
cat("Fig 1 : fig1_histogrammes.png\n")

# G2 : Boxplots par groupe ADA
p_box_list <- lapply(seq_along(vars_plot), function(i) {
  v <- vars_plot[i]
  ggplot(df_A, aes(x=ada_class, y=.data[[v]], fill=ada_class)) +
    geom_boxplot(outlier.colour="#FF4444", outlier.size=1.5, alpha=0.8) +
    scale_fill_manual(values=c("Normoglycémique"="#4ADE8088",
                               "Pré-diabétique"="#FFD70088",
                               "Diabétique"="#FF444488")) +
    labs(title=labels_plot[i], x="Groupe ADA", y=labels_plot[i]) +
    theme_minimal(base_size=10) +
    theme(legend.position="none",
          axis.text.x=element_text(angle=20, hjust=1),
          plot.title=element_text(hjust=0.5, face="bold"))
})
p_box <- ggarrange(plotlist=p_box_list, ncol=3, nrow=2)
ggsave("fig2_boxplots_ADA.png", p_box, width=14, height=9, dpi=300)
cat("Fig 2 : fig2_boxplots_ADA.png\n")

# G3 : Barres distributions
p_bar_ada  <- ggplot(df, aes(x=ada_class, fill=ada_class)) +
  geom_bar() +
  geom_text(stat="count",
            aes(label=paste0(after_stat(count),
                             "\n(",round(after_stat(count)/nrow(df)*100,1),"%)")),
            vjust=-0.3, size=3.2) +
  scale_fill_manual(values=c("Normoglycémique"="#4ADE80",
                             "Pré-diabétique"="#FFD700",
                             "Diabétique"="#FF4444",
                             "Hypoglycémie"="#38BDF8")) +
  labs(title="Classification ADA (n=462)",
       x="Groupe glycémique", y="Effectif") +
  theme_minimal(base_size=12) +
  theme(legend.position="none",
        axis.text.x=element_text(angle=15, hjust=1))

p_bar_cubt <- ggplot(df, aes(x=cubt_class, fill=cubt_class)) +
  geom_bar() +
  geom_text(stat="count",
            aes(label=paste0(after_stat(count),
                             "\n(",round(after_stat(count)/nrow(df)*100,1),"%)")),
            vjust=-0.3, size=3.5) +
  scale_fill_manual(values=c("Dans la plage"="#4ADE80","Hors plage"="#FF4444")) +
  labs(title="Classification CUBT (n=462)",
       x="Plage de référence", y="Effectif") +
  theme_minimal(base_size=12) + theme(legend.position="none")

ggsave("fig3_barres.png", ggarrange(p_bar_ada,p_bar_cubt,ncol=2),
       width=12, height=6, dpi=300)
cat("Fig 3 : fig3_barres.png\n")

# G4 : Boxplots indices dérivés
vars_idx    <- c("castelli","tg_hdl","ldl_hdl")
labels_idx  <- c("Castelli (CT/HDL)","TG/HDL","LDL/HDL")
p_idx_list  <- lapply(seq_along(vars_idx), function(i) {
  v <- vars_idx[i]
  ggplot(df_A, aes(x=ada_class, y=.data[[v]], fill=ada_class)) +
    geom_boxplot(outlier.colour="#FF4444", outlier.size=1.5, alpha=0.8) +
    scale_fill_manual(values=c("Normoglycémique"="#4ADE8088",
                               "Pré-diabétique"="#FFD70088",
                               "Diabétique"="#FF444488")) +
    labs(title=labels_idx[i], x="Groupe ADA", y=labels_idx[i]) +
    theme_minimal(base_size=10) +
    theme(legend.position="none",
          axis.text.x=element_text(angle=20, hjust=1),
          plot.title=element_text(hjust=0.5, face="bold"))
})
ggsave("fig4_boxplots_indices.png",
       ggarrange(plotlist=p_idx_list, ncol=3),
       width=14, height=5, dpi=300)
cat("Fig 4 : fig4_boxplots_indices.png\n")


# ============================================================
# ÉTAPE 3 — TESTS DE NORMALITÉ
# ============================================================
cat("\n=== TESTS DE NORMALITÉ ===\n")
cat("Effectifs par groupe ADA :\n")
print(table(df_A$ada_class))

test_norm <- function(data, var, groupe_var, groupe_levels) {
  bind_rows(lapply(groupe_levels, function(g) {
    x <- data[[var]][data[[groupe_var]] == g]; x <- x[!is.na(x)]
    n <- length(x)
    if (n < 50) { test <- shapiro.test(x);          methode <- "Shapiro-Wilk"  }
    else        { test <- nortest::lillie.test(x);   methode <- "KS-Lilliefors" }
    tibble(Variable=var, Groupe=g, n=n, Methode=methode,
           Stat=round(test$statistic,4), p=round(test$p.value,4),
           Normale=ifelse(test$p.value>0.05,"OUI","NON"))
  }))
}

norm_results_A <- bind_rows(lapply(vars_corr, function(v)
  test_norm(df_A, v, "ada_class",
            c("Normoglycémique","Pré-diabétique","Diabétique"))))
cat("\nNormalité Section A :\n")
print(as.data.frame(norm_results_A))

norm_global <- bind_rows(lapply(vars_corr, function(v) {
  x <- df[[v]]; x <- x[!is.na(x)]
  test <- nortest::lillie.test(x)
  tibble(Variable=v, n=length(x), Methode="KS-Lilliefors",
         Stat=round(test$statistic,4), p=round(test$p.value,4),
         Normale=ifelse(test$p.value>0.05,"OUI","NON"))
}))
cat("\nNormalité population globale :\n")
print(as.data.frame(norm_global))

cat("\nDécision : distribution non normale généralisée\n")
cat("→ Spearman (corrélations) · KW + Dunn (inter-groupes) · MW (Section B)\n")


# ============================================================
# ═══════════════════════════════════════════════════════════
# SECTION A — CLASSIFICATION ADA (3 groupes, n=460)
# ═══════════════════════════════════════════════════════════
# ============================================================

# ============================================================
# ÉTAPE 4A — CORRÉLATIONS HbA1c × PARAMÈTRES LIPIDIQUES
# ============================================================
cat("\n\n=== SECTION A — ÉTAPE 4A : CORRÉLATIONS HbA1c × LIPIDES ===\n")

corr_avec_IC <- function(x, y, methode="spearman") {
  if (methode=="pearson") {
    t <- cor.test(x, y, method="pearson", alternative="two.sided")
    tibble(r_rho=round(t$estimate,4), IC_inf=round(t$conf.int[1],4),
           IC_sup=round(t$conf.int[2],4), p=round(t$p.value,4),
           Methode="Pearson", Signif=ifelse(t$p.value<0.05,"*","ns"))
  } else {
    rho <- DescTools::SpearmanRho(x, y, conf.level=0.95)
    pv  <- cor.test(x, y, method="spearman", alternative="two.sided",
                    exact=FALSE)$p.value
    tibble(r_rho=round(rho[1],4), IC_inf=round(rho[2],4),
           IC_sup=round(rho[3],4), p=round(pv,4),
           Methode="Spearman", Signif=ifelse(pv<0.05,"*","ns"))
  }
}

# Matrice globale HbA1c × lipides
norm_hba1c <- nortest::lillie.test(df$hba1c)
cat("Normalité HbA1c (global) : p =", round(norm_hba1c$p.value,4),
    "→", ifelse(norm_hba1c$p.value>0.05,"Normale","Non normale"),"\n")

corr_matrix <- bind_rows(lapply(seq_along(vars_corr), function(i) {
  norm_y  <- nortest::lillie.test(df[[vars_corr[i]]])
  methode <- ifelse(norm_hba1c$p.value>0.05 & norm_y$p.value>0.05,
                    "pearson","spearman")
  res <- corr_avec_IC(df$hba1c, df[[vars_corr[i]]], methode)
  bind_cols(tibble(Variable=labels_corr[i]), res)
}))
cat("\nMatrice corrélations HbA1c × lipides (n=462) :\n")
print(as.data.frame(corr_matrix))

# Visualisation matrice inter-lipides
mat_corr <- cor(df[, vars_corr], method="spearman", use="complete.obs")
png("fig5_matrice_corr.png", width=1600, height=1400, res=200)
corrplot(mat_corr, method="color", type="upper", addCoef.col="white",
         number.cex=0.7, tl.col="black", tl.srt=45,
         col=colorRampPalette(c("#FF4444","#E0E0E0","#4ADE80"))(200),
         title="Matrice de corrélations (Spearman) — Paramètres lipidiques",
         mar=c(0,0,2,0))
dev.off()
cat("Fig 5 : fig5_matrice_corr.png\n")

# Corrélations stratifiées par groupe ADA (TOUS les lipides)
cat("\nCorrélations stratifiées HbA1c × lipides par groupe ADA :\n")
cat(sprintf("%-10s | %-20s : %8s %14s %8s %4s\n",
            "Variable","Groupe","rho","IC 95%","p","Signif"))
cat(strrep("-",75),"\n")
for (i in seq_along(vars_corr)) {
  for (g in c("Normoglycémique","Pré-diabétique","Diabétique")) {
    sub <- df_A %>% filter(ada_class == g)
    if (nrow(sub) < 5) next
    rho <- tryCatch(DescTools::SpearmanRho(sub$hba1c, sub[[vars_corr[i]]],
                                           conf.level=0.95), error=function(e) c(NA,NA,NA))
    pv  <- tryCatch(cor.test(sub$hba1c, sub[[vars_corr[i]]],
                             method="spearman", exact=FALSE)$p.value,
                    error=function(e) NA)
    cat(sprintf("%-10s | %-20s : rho=%6.3f [%6.3f;%6.3f] p=%.4f %s\n",
                labels_corr[i], g,
                rho[1], rho[2], rho[3], pv,
                ifelse(!is.na(pv) & pv<0.05,"*","ns")))
  }
}


# ==============================================================================
# CORRECTION 1 : FONCTION DE COMPARAISON (Étape 5A)
# ==============================================================================
compare_groupes <- function(data, var, groupe_var, norm_pval) {
  x_list <- split(data[[var]], data[[groupe_var]])
  
  if (norm_pval > 0.05) {
    # ANOVA de Welch
    model <- oneway.test(as.formula(paste(var, "~", groupe_var)), data=data, var.equal=FALSE)
    posthoc <- pairwise.t.test(data[[var]], data[[groupe_var]], p.adjust.method="bonferroni")
    methode <- "ANOVA Welch + post-hoc Bonferroni"
    stat <- model$statistic
    pv <- model$p.value
    
    cat(sprintf("\n%-12s : %s\n  Stat=%.4f  p=%.4f  %s\n", 
                var, methode, stat, pv, ifelse(pv<0.05,"DIFFÉRENCE SIGNIFICATIVE","ns")))
    print(posthoc$p.value)
    
  } else {
    # Kruskal-Wallis
    formule <- as.formula(paste(var, "~", groupe_var))
    kt <- kruskal.test(formule, data=data)
    
    # SOLUTION À L'ERREUR seq_len(nchar(fcLHS)) :
    # Utilisation de rstatix::dunn_test au lieu de FSA::dunnTest
    posthoc <- rstatix::dunn_test(data = data, formula = formule, p.adjust.method = "bonferroni")
    
    methode <- "Kruskal-Wallis + post-hoc Dunn (Bonferroni)"
    stat <- kt$statistic
    pv <- kt$p.value
    
    cat(sprintf("\n%-12s : %s\n  Stat=%.4f  p=%.4f  %s\n", 
                var, methode, stat, pv, ifelse(pv<0.05,"DIFFÉRENCE SIGNIFICATIVE","ns")))
    # Affichage allégé des résultats du post-hoc
    print(as.data.frame(posthoc[, c("group1", "group2", "statistic", "p", "p.adj.signif")]))
  }
  
  invisible(list(methode=methode, stat=stat, pvalue=pv, posthoc=posthoc))
}

# ==============================================================================
# CORRECTION 2 : APPLICATION ET NOMMAGE (Étape 5A suite)
# ==============================================================================
cat("\nComparaison des paramètres lipidiques entre les 3 groupes ADA :\n")

vars_anal <- c("hba1c","ct","ldl","hdl","tg","castelli","tg_hdl","ldl_hdl")

# Création de la liste avec la fonction sécurisée
comp_results <- lapply(vars_anal, function(v) {
  # Simulation du p global (à adapter selon votre objet norm_results_A existant)
  # Remplacer 0 par le vrai p_value si vous avez gardé la table norm_results_A
  norm_p <- 0 
  compare_groupes(df_A, v, "ada_class", norm_p)
})

# SOLUTION À L'ERREUR "Objet introuvable" :
# On s'assure d'assigner les noms avec le même vecteur qui a servi à créer la liste (vars_anal)
names(comp_results) <- vars_anal 

cat("\nTests inter-groupes réussis et enregistrés dans l'objet 'comp_results'.\n")

# G6 : Boxplot HbA1c avec test global
p_posthoc <- ggplot(df_A, aes(x=ada_ord, y=hba1c, fill=ada_ord)) +
  geom_boxplot(alpha=0.8) +
  geom_jitter(width=0.15, alpha=0.25, size=0.8, color="#C8D8E8") +
  scale_fill_manual(values=c("Normoglycémique"="#4ADE8099",
                             "Pré-diabétique"="#FFD70099",
                             "Diabétique"="#FF444499")) +
  stat_compare_means(method="kruskal.test",
                     label.y=max(df_A$hba1c, na.rm=TRUE)*1.05) +
  labs(title="HbA1c selon le groupe ADA (n=460)",
       x="Groupe glycémique", y="HbA1c (%)") +
  theme_minimal(base_size=12) + theme(legend.position="none")
ggsave("fig6_comparison_ADA.png", p_posthoc, width=9, height=6, dpi=300)
cat("\nFig 6 : fig6_comparison_ADA.png\n")


# ============================================================
# ÉTAPE 6A — COMPARAISON DES PROPORTIONS
# ============================================================
cat("\n=== ÉTAPE 6A : COMPARAISON DES PROPORTIONS ===\n")

test_proportion <- function(var_bin, groupe, data, label="") {
  tab      <- table(data[[var_bin]], data[[groupe]])
  cat(sprintf("\n%s × %s :\n", label, groupe)); print(tab)
  # [FIX] suppressWarnings pour le test intermédiaire d'effectifs
  chi_test <- suppressWarnings(chisq.test(tab))
  min_exp  <- min(chi_test$expected)
  if (min_exp >= 5) {
    cat(sprintf("Chi-2 : X2=%.4f  df=%d  p=%.4f  %s\n",
                chi_test$statistic, chi_test$parameter,
                chi_test$p.value, ifelse(chi_test$p.value<0.05,"*","ns")))
  } else {
    fisher <- fisher.test(tab, simulate.p.value=TRUE, B=10000)
    cat(sprintf("Fisher (effectif min=%.2f) : p=%.4f  %s\n",
                min_exp, fisher$p.value, ifelse(fisher$p.value<0.05,"*","ns")))
  }
}

vars_dich   <- c("ct_eleve","ct_bas","ldl_eleve","ldl_bas","hdl_eleve",
                 "hdl_bas","tg_eleve","tg_bas","castelli_patho",
                 "tg_hdl_patho","ldl_hdl_patho")
labels_dich <- c("CT élevé","CT bas","LDL élevé","LDL bas",
                 "HDL élevé","HDL bas","TG élevé","TG bas",
                 "Castelli patho.","TG/HDL pathol.","LDL/HDL pathol.")

for (i in seq_along(vars_dich))
  test_proportion(vars_dich[i], "ada_class", df_A, labels_dich[i])
test_proportion("sexe", "ada_class", df_A, "Sexe")


# ============================================================
# ÉTAPE 7A — MODÉLISATION 1 : RÉGRESSION LOGISTIQUE ORDINALE
# ============================================================
cat("\n\n=== ÉTAPE 7A — MODÉLISATION 1 : RÉGRESSION ORDINALE ===\n")

cat("
Note multicolinéarité [FIX-3] :
  CT = LDL + HDL + VLDL (relation algébrique partielle).
  Inclure CT et LDL ensemble génère VIF > 30 (invalide).
  → Modèle a1 : CT (sans LDL) | Modèle a2 : LDL (sans CT) | Modèle b : indices dérivés
")

# Mise à l'échelle + winsorisation TG
df_A_sc <- df_A %>%
  mutate(across(c(ct,ldl,hdl,tg,tg_wins,castelli,tg_hdl,tg_hdl_wins,
                  ldl_hdl,age),
                ~ scale(.x)[,1], .names="{.col}_z"))

# ─── Modèle 1a1 : CT + HDL + TG (sans LDL) ─────────────────
cat("\n--- Modèle 1a1 : CT + HDL + TG (sans LDL) ---\n")
m1a1 <- MASS::polr(
  ada_ord ~ ct_z + hdl_z + tg_wins_z + age_z + sexe,
  data=df_A_sc, Hess=TRUE, method="logistic")
coef_1a1 <- coef(summary(m1a1))
p_1a1    <- pnorm(abs(coef_1a1[,"t value"]), lower.tail=FALSE)*2
cat("\nCoefficients + p-values Modèle 1a1 :\n")
print(round(cbind(coef_1a1, "p"=p_1a1), 4))
ci_1a1  <- confint(m1a1)
OR_1a1  <- exp(cbind("OR"=coef(m1a1), ci_1a1))
cat("\nOdds Ratios — Modèle 1a1 :\n"); print(round(OR_1a1,4))
# [FIX-1] VIF direct sur polr
cat("\nVIF Modèle 1a1 :\n"); print(car::vif(m1a1))
cat("Brant test Modèle 1a1 :\n"); brant::brant(m1a1)

# ─── Modèle 1a2 : LDL + HDL + TG (sans CT) ─────────────────
cat("\n--- Modèle 1a2 : LDL + HDL + TG (sans CT) ---\n")
m1a2 <- MASS::polr(
  ada_ord ~ ldl_z + hdl_z + tg_wins_z + age_z + sexe,
  data=df_A_sc, Hess=TRUE, method="logistic")
coef_1a2 <- coef(summary(m1a2))
p_1a2    <- pnorm(abs(coef_1a2[,"t value"]), lower.tail=FALSE)*2
cat("\nCoefficients + p-values Modèle 1a2 :\n")
print(round(cbind(coef_1a2, "p"=p_1a2), 4))
ci_1a2  <- confint(m1a2)
OR_1a2  <- exp(cbind("OR"=coef(m1a2), ci_1a2))
cat("\nOdds Ratios — Modèle 1a2 :\n"); print(round(OR_1a2,4))
cat("\nVIF Modèle 1a2 :\n"); print(car::vif(m1a2))
cat("Brant test Modèle 1a2 :\n"); brant::brant(m1a2)

# ─── Modèle 1b : indices dérivés ────────────────────────────
cat("\n--- Modèle 1b : indices dérivés ---\n")
m1b <- MASS::polr(
  ada_ord ~ castelli_z + tg_hdl_wins_z + ldl_hdl_z + age_z + sexe,
  data=df_A_sc, Hess=TRUE, method="logistic")
coef_1b <- coef(summary(m1b))
p_1b    <- pnorm(abs(coef_1b[,"t value"]), lower.tail=FALSE)*2
cat("\nCoefficients + p-values Modèle 1b :\n")
print(round(cbind(coef_1b, "p"=p_1b), 4))
ci_1b   <- confint(m1b)
OR_1b   <- exp(cbind("OR"=coef(m1b), ci_1b))
cat("\nOdds Ratios — Modèle 1b :\n"); print(round(OR_1b,4))
cat("\nVIF Modèle 1b :\n"); print(car::vif(m1b))
cat("Brant test Modèle 1b :\n"); brant::brant(m1b)

cat("\nCritères d'information — Modèles 1 :\n")
cat("AIC 1a1:", AIC(m1a1),"| BIC:", BIC(m1a1),"\n")
cat("AIC 1a2:", AIC(m1a2),"| BIC:", BIC(m1a2),"\n")
cat("AIC 1b :", AIC(m1b), "| BIC:", BIC(m1b),  "\n")

# ─── Modèle alternatif multinomial (si Brant violé) ─────────
cat("\n--- Modèle alternatif : Régression multinomiale (nnet) ---\n")
cat("(Robuste à la violation de l'hypothèse de proportionnalité)\n")
m1_multi <- nnet::multinom(
  ada_ord ~ ct_z + hdl_z + tg_wins_z + age_z + sexe,
  data=df_A_sc, trace=FALSE)
cat("\nCoefficients multinomiaux :\n"); print(summary(m1_multi))
# Calcul p-values
z_multi <- summary(m1_multi)$coefficients / summary(m1_multi)$standard.errors
p_multi  <- (1 - pnorm(abs(z_multi))) * 2
cat("\np-values multinomiales :\n"); print(round(p_multi,4))
OR_multi <- exp(coef(m1_multi))
cat("\nOR multinomiaux :\n"); print(round(OR_multi,4))
cat("AIC multinomial:", AIC(m1_multi), "\n")


# ============================================================
# ÉTAPE 8A — MODÉLISATION 2 : RÉGRESSION BINAIRE (ADA)
# ============================================================
cat("\n\n=== ÉTAPE 8A — MODÉLISATION 2 : RÉGRESSION BINAIRE (ADA) ===\n")
cat("Non pathologique (Normo) vs Pathologique (Prédiab + Diab)\n")
print(table(df_A_bin$ada_bin))

df_A_bin_sc <- df_A_bin %>%
  mutate(across(c(ct,ldl,hdl,tg,tg_wins,castelli,tg_hdl,tg_hdl_wins,
                  ldl_hdl,age),
                ~ scale(.x)[,1], .names="{.col}_z"),
         ada_bin_num = as.numeric(ada_bin)-1)

# [FIX-7] nagelkerke_r2 corrigée
nagelkerke_r2 <- function(model) {
  null_ll <- as.numeric(logLik(update(model,.~1)))
  full_ll <- as.numeric(logLik(model))
  n       <- nobs(model)
  r2_cs   <- 1 - exp((2/n)*(null_ll-full_ll))
  r2_max  <- 1 - exp((2/n)*null_ll)
  round(r2_cs/r2_max, 4)
}

# ─── Modèle 2a1 : CT + HDL + TG [FIX-3] ────────────────────
cat("\n--- Modèle 2a1 : CT + HDL + TG (sans LDL) ---\n")
m2a1 <- glm(ada_bin ~ ct_z + hdl_z + tg_wins_z + age_z + sexe,
            data=df_A_bin_sc, family=binomial(link="logit"))
print(summary(m2a1))
OR_2a1 <- exp(cbind("OR"=coef(m2a1), suppressMessages(confint(m2a1))))
cat("\nOR — Modèle 2a1 :\n"); print(round(OR_2a1,4))
cat("\nVIF 2a1 :\n"); print(car::vif(m2a1))

# ─── Modèle 2a2 : LDL + HDL + TG ───────────────────────────
cat("\n--- Modèle 2a2 : LDL + HDL + TG (sans CT) ---\n")
m2a2 <- glm(ada_bin ~ ldl_z + hdl_z + tg_wins_z + age_z + sexe,
            data=df_A_bin_sc, family=binomial(link="logit"))
print(summary(m2a2))
OR_2a2 <- exp(cbind("OR"=coef(m2a2), suppressMessages(confint(m2a2))))
cat("\nOR — Modèle 2a2 :\n"); print(round(OR_2a2,4))
cat("\nVIF 2a2 :\n"); print(car::vif(m2a2))

# ─── Modèle 2b : indices dérivés ────────────────────────────
cat("\n--- Modèle 2b : indices dérivés ---\n")
m2b <- glm(ada_bin ~ castelli_z + tg_hdl_wins_z + ldl_hdl_z + age_z + sexe,
           data=df_A_bin_sc, family=binomial(link="logit"))
print(summary(m2b))
OR_2b <- exp(cbind("OR"=coef(m2b), suppressMessages(confint(m2b))))
cat("\nOR — Modèle 2b :\n"); print(round(OR_2b,4))
cat("\nVIF 2b :\n"); print(car::vif(m2b))

# ─── Modèle 2c : Firth (robuste séparation) ─────────────────
cat("\n--- Modèle 2c : Firth (alternative robuste — séparation parfaite) ---\n")
m2c <- logistf(ada_bin ~ ct_z + hdl_z + tg_wins_z + age_z + sexe,
               data=df_A_bin_sc, firth=TRUE)
cat("\nRésumé Firth 2c :\n"); print(summary(m2c))

# ─── Diagnostics modèle 2b (meilleur AIC) ───────────────────
cat("\n--- Diagnostics — Modèle 2b ---\n")
hl_2b <- ResourceSelection::hoslem.test(
  df_A_bin_sc$ada_bin_num, fitted(m2b), g=10)
cat(sprintf("Hosmer-Lemeshow : X2=%.4f  p=%.4f  %s\n",
            hl_2b$statistic, hl_2b$p.value,
            ifelse(hl_2b$p.value>0.05,"Bon ajustement","Ajustement insuffisant")))

roc_2b <- pROC::roc(df_A_bin_sc$ada_bin_num, fitted(m2b), quiet=TRUE)
ci_roc_2b <- pROC::ci(roc_2b)
cat(sprintf("AUC Modèle 2b : %.4f [%.4f;%.4f]\n",
            as.numeric(pROC::auc(roc_2b)), ci_roc_2b[1], ci_roc_2b[3]))

png("fig7_ROC_modele2b.png", width=900, height=900, res=150)
plot(roc_2b, col="#12B3C4", lwd=2,
     main=paste0("Courbe ROC — Modèle 2b ADA (AUC = ",
                 round(pROC::auc(roc_2b),3),")"),
     print.auc=TRUE)
dev.off()
cat("Fig 7 : fig7_ROC_modele2b.png\n")

# Matrice de confusion
seuil_opt <- pROC::coords(roc_2b, "best", ret="threshold")$threshold
pred_class <- ifelse(fitted(m2b) >= seuil_opt, "Pathologique","Non pathologique")
cat("\nMatrice de confusion (seuil optimal =", round(seuil_opt,3),") :\n")
print(table(Pred=pred_class, Obs=df_A_bin_sc$ada_bin))

cat("Pseudo-R2 Nagelkerke — 2a1:", nagelkerke_r2(m2a1), "\n")
cat("Pseudo-R2 Nagelkerke — 2a2:", nagelkerke_r2(m2a2), "\n")
cat("Pseudo-R2 Nagelkerke — 2b :", nagelkerke_r2(m2b), "\n")
cat("AIC 2a1:", AIC(m2a1),"| 2a2:", AIC(m2a2),"| 2b:", AIC(m2b),"\n")


# ============================================================
# ═══════════════════════════════════════════════════════════
# SECTION B — CLASSIFICATION CUBT (2 groupes, n=462)
# ═══════════════════════════════════════════════════════════
# ============================================================

cat("\n\n=== SECTION B — CLASSIFICATION CUBT ===\n")
cat("Distribution :\n"); print(table(df_B$cubt_class))


# ============================================================
# ÉTAPE 4B — COMPARAISON INTER-GROUPES CUBT
# ============================================================
cat("\n=== ÉTAPE 4B : COMPARAISON INTER-GROUPES CUBT ===\n")
cat("[FIX-2] Utilise vars_corr (7 lipides, hba1c exclus)\n\n")

compare_2groupes <- function(data, var, groupe_var, norm_p, label="") {
  g_levels <- levels(data[[groupe_var]])
  x1 <- data[[var]][data[[groupe_var]]==g_levels[1]]; x1 <- x1[!is.na(x1)]
  x2 <- data[[var]][data[[groupe_var]]==g_levels[2]]; x2 <- x2[!is.na(x2)]
  cat(sprintf("\n%s :\n", label))
  cat(sprintf("  %s : n=%d | med=%.3f | moy=%.3f\n",
              g_levels[1], length(x1), median(x1), mean(x1)))
  cat(sprintf("  %s  : n=%d | med=%.3f | moy=%.3f\n",
              g_levels[2], length(x2), median(x2), mean(x2)))
  if (norm_p > 0.05) {
    test <- t.test(x1, x2, alternative="two.sided",
                   var.equal=FALSE, conf.level=0.95)
    cat(sprintf("  t-test Welch : t=%.4f  df=%.1f  p=%.4f  IC95%%=[%.4f;%.4f]  %s\n",
                test$statistic, test$parameter, test$p.value,
                test$conf.int[1], test$conf.int[2],
                ifelse(test$p.value<0.05,"*","ns")))
  } else {
    test <- wilcox.test(x1, x2, alternative="two.sided",
                        conf.int=TRUE, conf.level=0.95, exact=FALSE)
    cat(sprintf("  Mann-Whitney : W=%.1f  p=%.4f  IC95%%=[%.4f;%.4f]  %s\n",
                test$statistic, test$p.value,
                test$conf.int[1], test$conf.int[2],
                ifelse(test$p.value<0.05,"*","ns")))
  }
}

# [FIX-2] Ajouter HbA1c manuellement en premier
cat("\nhba1c :")
cat("\n  Dans la plage : n=",sum(df_B$cubt_class=="Dans la plage"),
    "| med=", median(df_B$hba1c[df_B$cubt_class=="Dans la plage"]),
    "| moy=", round(mean(df_B$hba1c[df_B$cubt_class=="Dans la plage"]),3))
cat("\n  Hors plage    : n=",sum(df_B$cubt_class=="Hors plage"),
    "| med=", median(df_B$hba1c[df_B$cubt_class=="Hors plage"]),
    "| moy=", round(mean(df_B$hba1c[df_B$cubt_class=="Hors plage"]),3),"\n")
wt_hba1c <- wilcox.test(hba1c~cubt_bin, data=df_B, exact=FALSE,
                        conf.int=TRUE, conf.level=0.95)
cat(sprintf("  Mann-Whitney : W=%.1f  p=%.4f  IC95%%=[%.4f;%.4f]  %s\n",
            wt_hba1c$statistic, wt_hba1c$p.value,
            wt_hba1c$conf.int[1], wt_hba1c$conf.int[2],
            ifelse(wt_hba1c$p.value<0.05,"*","ns")))

# [FIX-2] Boucle correcte sur vars_corr et labels_corr
for (i in seq_along(vars_corr)) {
  norm_p_b <- nortest::lillie.test(df_B[[vars_corr[i]]])$p.value
  compare_2groupes(df_B, vars_corr[i], "cubt_bin",
                   norm_p_b, labels_corr[i])
}

# G8 : Boxplots Section B
p_box_B <- lapply(seq_along(vars_plot), function(i) {
  v <- vars_plot[i]
  ggplot(df_B, aes(x=cubt_class, y=.data[[v]], fill=cubt_class)) +
    geom_boxplot(outlier.colour="#FF4444", outlier.size=1.5, alpha=0.8) +
    scale_fill_manual(values=c("Dans la plage"="#4ADE8088",
                               "Hors plage"="#FF444488")) +
    labs(title=labels_plot[i], x="Plage CUBT", y=labels_plot[i]) +
    theme_minimal(base_size=10) +
    theme(legend.position="none",
          plot.title=element_text(hjust=0.5, face="bold"))
})
ggsave("fig8_boxplots_CUBT.png",
       ggarrange(plotlist=p_box_B, ncol=3, nrow=2),
       width=14, height=9, dpi=300)
cat("\nFig 8 : fig8_boxplots_CUBT.png\n")


# ============================================================
# ÉTAPE 5B — COMPARAISON DES PROPORTIONS (CUBT)
# ============================================================
cat("\n=== ÉTAPE 5B : COMPARAISON DES PROPORTIONS — CUBT ===\n")
for (i in seq_along(vars_dich))
  test_proportion(vars_dich[i], "cubt_class", df_B, labels_dich[i])
test_proportion("sexe","cubt_class", df_B, "Sexe")


# ============================================================
# ÉTAPE 6B — MODÉLISATION 3 : RÉGRESSION BINAIRE (CUBT)
# ============================================================
cat("\n\n=== ÉTAPE 6B — MODÉLISATION 3 : RÉGRESSION BINAIRE (CUBT) ===\n")
cat("Note : HbA1c est la variable de définition du groupe CUBT\n")
cat("→ HbA1c EXCLU des prédicteurs (tautologie)\n\n")
print(table(df_B$cubt_bin))

df_B_sc <- df_B %>%
  mutate(across(c(ct,ldl,hdl,tg,tg_wins,castelli,tg_hdl,tg_hdl_wins,
                  ldl_hdl,age),
                ~ scale(.x)[,1], .names="{.col}_z"),
         cubt_num = as.numeric(cubt_bin)-1)

# ─── Modèle 3a1 : CT + HDL + TG ────────────────────────────
cat("\n--- Modèle 3a1 : CT + HDL + TG (sans LDL) ---\n")
m3a1 <- glm(cubt_bin ~ ct_z + hdl_z + tg_wins_z + age_z + sexe,
            data=df_B_sc, family=binomial(link="logit"))
print(summary(m3a1))
OR_3a1 <- exp(cbind("OR"=coef(m3a1), suppressMessages(confint(m3a1))))
cat("\nOR — Modèle 3a1 :\n"); print(round(OR_3a1,4))
cat("VIF 3a1 :\n"); print(car::vif(m3a1))

# ─── Modèle 3a2 : LDL + HDL + TG ───────────────────────────
cat("\n--- Modèle 3a2 : LDL + HDL + TG (sans CT) ---\n")
m3a2 <- glm(cubt_bin ~ ldl_z + hdl_z + tg_wins_z + age_z + sexe,
            data=df_B_sc, family=binomial(link="logit"))
print(summary(m3a2))
OR_3a2 <- exp(cbind("OR"=coef(m3a2), suppressMessages(confint(m3a2))))
cat("\nOR — Modèle 3a2 :\n"); print(round(OR_3a2,4))

# ─── Modèle 3b : indices dérivés ────────────────────────────
cat("\n--- Modèle 3b : indices dérivés ---\n")
m3b <- glm(cubt_bin ~ castelli_z + tg_hdl_wins_z + ldl_hdl_z + age_z + sexe,
           data=df_B_sc, family=binomial(link="logit"))
print(summary(m3b))
OR_3b <- exp(cbind("OR"=coef(m3b), suppressMessages(confint(m3b))))
cat("\nOR — Modèle 3b :\n"); print(round(OR_3b,4))

# ─── Diagnostics modèle 3a1 ─────────────────────────────────
cat("\n--- Diagnostics — Modèle 3a1 ---\n")
hl_3a1 <- ResourceSelection::hoslem.test(df_B_sc$cubt_num, fitted(m3a1), g=10)
cat(sprintf("Hosmer-Lemeshow : X2=%.4f  p=%.4f  %s\n",
            hl_3a1$statistic, hl_3a1$p.value,
            ifelse(hl_3a1$p.value>0.05,"Bon ajustement","Ajustement insuffisant")))

roc_3a1 <- pROC::roc(df_B_sc$cubt_num, fitted(m3a1), quiet=TRUE)
ci_roc_3  <- pROC::ci(roc_3a1)
cat(sprintf("AUC Modèle 3a1 : %.4f [%.4f;%.4f]\n",
            as.numeric(pROC::auc(roc_3a1)), ci_roc_3[1], ci_roc_3[3]))

png("fig9_ROC_modele3a1.png", width=900, height=900, res=150)
plot(roc_3a1, col="#F97316", lwd=2,
     main=paste0("Courbe ROC — Modèle 3a1 CUBT (AUC = ",
                 round(pROC::auc(roc_3a1),3),")"),
     print.auc=TRUE)
dev.off()
cat("Fig 9 : fig9_ROC_modele3a1.png\n")

cat("Pseudo-R2 Nagelkerke — 3a1:", nagelkerke_r2(m3a1), "\n")
cat("Pseudo-R2 Nagelkerke — 3a2:", nagelkerke_r2(m3a2), "\n")
cat("Pseudo-R2 Nagelkerke — 3b :", nagelkerke_r2(m3b), "\n")
cat("AIC 3a1:", AIC(m3a1),"| 3a2:", AIC(m3a2),"| 3b:", AIC(m3b),"\n")


# ============================================================
# RÉCAPITULATIF FINAL
# ============================================================
cat("\n\n========================================================\n")
cat("RÉCAPITULATIF — TOUS LES MODÈLES (version corrigée)\n")
cat("========================================================\n")

recap <- tribble(
  ~Modele, ~Type,         ~Section, ~Predicteurs,                 ~AIC,
  "M1a1", "Ordinal",      "A",      "CT + HDL + TG (wins)",       AIC(m1a1),
  "M1a2", "Ordinal",      "A",      "LDL + HDL + TG (wins)",      AIC(m1a2),
  "M1b",  "Ordinal",      "A",      "Indices dérivés",             AIC(m1b),
  "M1m",  "Multinomial",  "A",      "CT + HDL + TG (wins)",       AIC(m1_multi),
  "M2a1", "Binaire",      "A",      "CT + HDL + TG (wins)",       AIC(m2a1),
  "M2a2", "Binaire",      "A",      "LDL + HDL + TG (wins)",      AIC(m2a2),
  "M2b",  "Binaire",      "A",      "Indices dérivés",             AIC(m2b),
  "M3a1", "Binaire",      "B",      "CT + HDL + TG (wins)",       AIC(m3a1),
  "M3a2", "Binaire",      "B",      "LDL + HDL + TG (wins)",      AIC(m3a2),
  "M3b",  "Binaire",      "B",      "Indices dérivés",             AIC(m3b)
)
print(as.data.frame(recap))
cat("\n========================================================\n")
cat("Analyse terminée — v2 corrigée\n")
cat("Figures produites : fig1 à fig9\n")
cat("========================================================\n")

