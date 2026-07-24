# =============================================================================
# Script 16: Sync result.docx — add Figure 1B cumulative incidence reference
# =============================================================================

suppressPackageStartupMessages({ library(officer) })

doc_path <- "analysis/result.docx"
doc      <- read_docx(doc_path)

TNR  <- "Times New Roman"
fp_bd <- fp_text(font.family = TNR, font.size = 11, bold = TRUE)
fp_it <- fp_text(font.family = TNR, font.size = 11, italic = TRUE)
fp_nm <- fp_text(font.family = TNR, font.size = 11)

# Locate Figure 1 line and insert Figure 1B immediately after
doc <- cursor_reach(doc, keyword = "Figure 2 — see analysis/figures/Figure2_KM_Equipoise")
doc <- cursor_forward(doc)

doc <- body_add_fpar(doc,
  fpar(
    ftext("Figure 1B.  ", fp_bd),
    ftext(paste0(
      "Cumulative Incidence of 28-Day All-Cause Death by Treatment Group ",
      "(Equipoise Subgroup, SOFA 8-14; n = 1,952). ",
      "The cumulative incidence curve (1 - Kaplan-Meier estimator) presents ",
      "the same IPCW-weighted survival data as Figure 1, reframed to show the ",
      "probability of death accumulating over time. The teal curve represents ",
      "the Early CRRT group and the vermillion curve the Deferred CRRT group, ",
      "with 95% confidence bands. The two representations are mathematically ",
      "equivalent; the cumulative incidence format is included to facilitate ",
      "direct visual comparison with the CIF curves for the secondary outcome ",
      "(Figure 3). IPCW-weighted Cox p-value is provided for reference."
    ), fp_it),
    fp_p = fp_par(padding.top = 4, padding.bottom = 6)
  ),
  style = "Normal"
)

# Also update Suppl. Fig S1 reference to include the corresponding full-cohort version
doc <- cursor_reach(doc, keyword = "Figure S1 — see analysis/figures/FigS1_KM_FullCohort")
doc <- cursor_forward(doc)

doc <- body_add_fpar(doc,
  fpar(
    ftext("[Supplementary Figure S1b — see analysis/figures/FigS1b_CumInc_Full.pdf]", fp_nm),
    fp_p = fp_par(padding.top = 0, padding.bottom = 6)
  ),
  style = "Normal"
)

# Update Supplementary Materials section — add S1b entry
doc <- cursor_reach(doc, keyword = "Supplementary Figure S1.")
doc <- cursor_forward(doc)

doc <- body_add_fpar(doc,
  fpar(
    ftext("Supplementary Figure S1b.  ",
          fp_text(font.family = TNR, font.size = 9, bold = TRUE)),
    ftext(paste0(
      "Cumulative incidence of 28-day all-cause death, full analytic cohort ",
      "(n = 3,112). Companion to Supplementary Figure S1; same IPCW-weighted ",
      "data presented as 1 - KM. ",
      "[See: analysis/figures/FigS1b_CumInc_Full.pdf]"
    ), fp_text(font.family = TNR, font.size = 9, color = "grey30")),
    fp_p = fp_par(padding.top = 3, padding.bottom = 3)
  ),
  style = "Normal"
)

print(doc, target = doc_path)
cat("result.docx updated with Figure 1B reference.\n")
