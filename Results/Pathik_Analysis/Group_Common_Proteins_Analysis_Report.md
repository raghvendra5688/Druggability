# Pathik Compound Group — Common Protein Target Analysis

**Source file:** `Group_Common_Proteins_All_Thresholds.csv`
**Method:** For each compound group, DrugBank analogs were identified by Tanimoto similarity (threshold ≥ 0.3, top-20 per compound). Proteins ranking in the top-k (by descending druggability score) for **all** common DrugBank analogs were intersected to yield group-level target predictions. Analysis was run at three thresholds: top-200, top-500, and top-1000.

---

## Group Overview

| Group | Compounds | # Compounds | # Common DrugBank Analogs |
|---|---|---|---|
| 1 | Compound_2, Compound_3, Compound_31, Compound_53 | 4 | 14 |
| 2 | Compound_5, Compound_6 | 2 | 10 |
| 3 | Compound_8, Compound_4, Compound_7 | 3 | 2 |
| 4 | Compound_13, Compound_19, Compound_28 | 3 | 6 |
| 5 | Novel1, Novel2, Novel3 | 3 | 10 |
| 6 | Novel4, Novel5 | 2 | 0 |

Group 6 produced no common analogs at the Tanimoto threshold and is excluded from downstream analysis.

---

## Common Protein Counts per Group

| Group | Top-200 | Top-500 | Top-1000 |
|---|---|---|---|
| 1 | 24 | 55 | 100 |
| 2 | 24 | 59 | 118 |
| 3 | 76 | 193 | 393 |
| 4 | 25 | 65 | 147 |
| 5 | 24 | 49 | 97 |
| 6 | 0 | 0 | 0 |

Group 3 consistently yields the largest intersection, likely because it has only 2 common DrugBank analogs (Tariquidar and Pagoclone), making the intersection less stringent.

---

## Known Drug Targets in Top Hits

### Proteins appearing at top-200 (most stringent — common to all groups 1–5)

These proteins are shared across nearly every group at the strictest threshold and represent the most robust predictions.

| Protein | Full Name | Known Drug Target? | Example Drugs |
|---|---|---|---|
| **AR** | Androgen Receptor | Yes | Enzalutamide, bicalutamide, apalutamide |
| **TNF** | Tumour Necrosis Factor | Yes | Adalimumab, etanercept, infliximab |
| **EGF** | Epidermal Growth Factor | Yes (via EGFR) | Cetuximab, erlotinib, gefitinib |
| **F2** | Prothrombin (Factor II) | Yes | Dabigatran, argatroban |
| **F3** | Tissue Factor (Factor III) | Yes | Factor VIIa inhibitors |
| **F5** | Coagulation Factor V | Yes | Anticoagulation pathway target |
| **F7** | Coagulation Factor VII | Yes | Eptacog alfa |
| **F8** | Coagulation Factor VIII | Yes | Emicizumab, recombinant FVIII |
| **F9** | Coagulation Factor IX | Yes | Fitusiran, recombinant FIX |
| **F10** | Coagulation Factor X | Yes | Rivaroxaban, apixaban, edoxaban |
| **F11** | Coagulation Factor XI | Yes | Abelacimab (in trials) |
| **C3** | Complement Component 3 | Yes | Pegcetacoplan (FDA approved) |
| **C5** | Complement Component 5 | Yes | Eculizumab, ravulizumab (FDA approved) |
| **C2, C6, C7, C9** | Complement pathway | Partial | C6/C7/C9 targeted experimentally |
| **NF1** | Neurofibromin 1 | Indirect | Loss activates RAS/MEK — targeted via MEK inhibitors |
| SI, CS, GC, HP, MB, ALB, TH, CDA, KL, TG | Various | Generally no | Housekeeping / plasma proteins — likely scoring noise |

---

### Additional known targets emerging at top-500

| Protein | Full Name | Known Drug Target? | Example Drugs / Context |
|---|---|---|---|
| **BCL2** | B-cell lymphoma 2 | Yes | Venetoclax, navitoclax |
| **MYC** | MYC proto-oncogene | Emerging | BET inhibitors (indirect); direct inhibitors in trials |
| **AKT1** | AKT serine/threonine kinase 1 | Yes | Ipatasertib, capivasertib |
| **MET** | MET receptor tyrosine kinase | Yes | Crizotinib, cabozantinib, tepotinib |
| **ALK** | Anaplastic lymphoma kinase | Yes | Crizotinib, alectinib, lorlatinib |
| **HRAS / NRAS** | RAS GTPases | Yes | Pathway targeted via MEK/ERK inhibitors |
| **JUN / FOS** | AP-1 transcription factors | Indirect | Downstream of many oncogenic pathways |
| **NFKBIA** | IκBα (NF-κB inhibitor) | Indirect | NF-κB pathway — bortezomib targets upstream |
| **RB1** | Retinoblastoma protein | Indirect (Group 3) | Loss drives CDK4/6 dependence — palbociclib context |

---

### Additional known targets emerging at top-1000

| Protein | Full Name | Known Drug Target? | Example Drugs / Context |
|---|---|---|---|
| **MTOR** | mTOR kinase | Yes | Rapamycin, everolimus, temsirolimus |
| **SRC** | SRC tyrosine kinase | Yes | Dasatinib, bosutinib |
| **MAPK1 / MAPK3** | ERK1/2 | Yes | Ulixertinib, MK-8353 |
| **STAT3** | Signal transducer and activator 3 | Emerging | Multiple STAT3 inhibitors in trials |
| **CTNNB1** | β-catenin | Emerging | Wnt pathway inhibitors in trials |
| **CDK1 / CDK5** | Cyclin-dependent kinases | Yes | CDK inhibitors (palbociclib, ribociclib class) |
| **TP53** | Tumour protein p53 | Indirect (Group 3) | MDM2 inhibitors restore TP53 activity |
| **IKBKG** | NF-κB essential modulator | Indirect (Group 3) | Upstream of NF-κB survival signalling |
| **PTEN** | Phosphatase and tensin homolog | Indirect | Loss activates PI3K/AKT; targeted via AKT inhibitors |
| **FGFR1** | FGF receptor 1 | Yes (Group 5) | Erdafitinib, pemigatinib |
| **ESR1** | Estrogen receptor 1 | Yes (Group 3) | Tamoxifen, fulvestrant, elacestrant |
| **PLK1** | Polo-like kinase 1 | Yes (Group 3) | Volasertib, onvansertib |
| **KRAS** | KRAS GTPase | Yes (Group 3) | Sotorasib, adagrasib |

---

## Key Findings by Group

### Group 1 (Compound_2, _3, _31, _53)
- Top-200 hits dominated by coagulation factors and complement proteins alongside AR and EGF.
- BCL2 and AKT1 emerge at top-1000, suggesting these compounds may engage apoptosis and survival signalling.
- MYC, JUN, FOS, TNF, REL also appear at top-1000, pointing to broad transcriptional/inflammatory pathway involvement.

### Group 2 (Compound_5, Compound_6)
- Similar coagulation/complement profile at top-200 as Group 1.
- At top-1000, BCL2, HRAS, MYC, MET, JUN, FOS, and MAX appear, indicating overlap with oncogenic signalling networks.
- NFKBIA at top-500/1000 suggests NF-κB pathway relevance.

### Group 3 (Compound_8, _4, _7)
- Largest intersection at all thresholds (76 / 193 / 393 proteins) due to only 2 DrugBank analogs.
- Most diverse target profile: BCL2 (top-500), RB1, CDK5, AKT1 (top-500), and TP53, IKBKG, CTNNB1, KRAS, NRAS, MAPK1/3, ESR1, PLK1, PTEN, STAT3 at top-1000.
- Suggests these compounds may act on broad oncogenic programmes spanning cell cycle, apoptosis, RAS/MAPK, and Wnt signalling.

### Group 4 (Compound_13, _19, _28)
- TNF prominent at top-200 alongside coagulation/complement targets and AR.
- ALK and MET appear at top-500, suggesting receptor tyrosine kinase relevance.
- BCL2 and MTOR emerge at top-1000, consistent with thalidomide/lenalidomide/pomalidomide analogs (immunomodulatory drugs with known MTOR and anti-apoptotic activity).

### Group 5 (Novel1, Novel2, Novel3)
- BCL2 appears early — at top-500 — making it the most prominent actionable target for this group.
- MTOR, SRC, MAPK3, MET, FGFR1, AKT1, HRAS also appear at top-1000.
- STAT3 appears at both top-500 and top-1000, consistent with vinca alkaloid-like DrugBank analogs that disrupt mitotic/survival signalling.

### Group 6 (Novel4, Novel5)
- No common DrugBank analogs found at Tanimoto threshold ≥ 0.3. No protein predictions possible. Consider lowering the threshold or expanding the reference database.

---

## Summary Assessment

The predicted top targets are **largely biologically plausible**:

- **Coagulation factors and complement proteins** dominate at top-200 across all groups — these are well-validated drug target families and their consistent appearance suggests the scoring captures known pharmacology.
- **AR and TNF** appear at top-200 in every group — both are major drug target classes, though their relevance may reflect broad scoring rather than group-specific biology.
- **BCL2** correctly emerges in all groups (top-500 or top-1000), consistent with the senolytic focus of the project and the known activity of venetoclax/navitoclax.
- **MTOR, AKT1, SRC, MAPK1/3, MET, FGFR1** — all clinically actionable kinases — appear at top-1000, supporting the druggability framework.
- **TP53, PTEN, RB1** are tumour suppressors and are not directly druggable but their appearance reflects oncogenic context.
- Housekeeping proteins (ALB, SI, GAPDH, CS, HP, MB) appearing consistently are likely artifacts of the scoring methodology rather than true target predictions.
