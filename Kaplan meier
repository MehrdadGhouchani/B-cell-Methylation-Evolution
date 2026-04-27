# Kaplan-Meier Survival Curves: TTT and OS
# =============================================
# Project: B cell lymphomas DNA methylation analysis
# April 2026

library(survival)
library(survminer)

# -----------------------------------------------------------------------------
# 1. Load and Prepare Clinical Data
# -----------------------------------------------------------------------------

CLL <- read.csv("clinic.csv")

# Parse dates
CLL$last_followup    <- as.Date(CLL$last_followup,    "%m/%d/%Y")
CLL$first_treatment  <- as.Date(CLL$first_treatment,  "%m/%d/%Y")
CLL$sampling         <- as.Date(CLL$sampling,         "%m/%d/%Y")

# Compute OS and TTT from sampling date (in years)
CLL$OS_sampling  <- as.numeric(CLL$last_followup   - CLL$sampling) / 365.25
CLL$TTT_sampling <- as.numeric(CLL$first_treatment - CLL$sampling) / 365.25

# Patients never treated: use OS as TTT (censored at last follow-up)
CLL$TTT_sampling[CLL$Need_for_treatment == 0] <- CLL$OS_sampling[CLL$Need_for_treatment == 0]

# -----------------------------------------------------------------------------
# 2. Filter Samples
# -----------------------------------------------------------------------------

# Remove pretreated samples (treatment started before sampling)
CLL <- CLL[CLL$TTT_sampling >= 0, ]

# Remove competing events: patients who died without needing treatment
# (applies to TTT only — do not apply this filter for OS analysis)
CLL <- CLL[!(CLL$Death == 1 & CLL$Need_for_treatment == 0), ]

# -----------------------------------------------------------------------------
# 3. Kaplan-Meier Curves
# -----------------------------------------------------------------------------

# 1. Overall Survival (OS)
surv_OS <- Surv(time = CLL$OS_sampling, event = CLL$Death)
fit_OS  <- survfit(surv_OS ~ CLL_epitype, data = CLL)

ggsurvplot(fit_OS,
           data           = CLL,
           risk.table     = TRUE,
           risk.table.col = "strata",
           pval           = TRUE,
           title          = "Overall Survival - CLL epitypes",
           ylab           = "Survival Probability")

# 2. Time to Treatment (TTT)
surv_TTT <- Surv(time = CLL$TTT_sampling, event = CLL$Need_for_treatment)
fit_TTT  <- survfit(surv_TTT ~ CLL_epitype, data = CLL)

ggsurvplot(fit_TTT,
           data           = CLL,
           risk.table     = TRUE,
           risk.table.col = "strata",
           pval           = TRUE,
           title          = "Time to Treatment - CLL epitypes",
           ylab           = "% Untreated")
