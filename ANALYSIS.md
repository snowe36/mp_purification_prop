# Analysis contract

**Status: pre-registered.** Written before Gate 0 census numbers and before any model is fit. Amendments are logged at the bottom.

The phenotype of interest is **HTS-studyable membrane protein material**, not fluorescence and not a published structure. Those are different measurements. They are never pooled.

---

## 1. Questions (never mixed)

### A — cellular expression

Did the construct accumulate in the membrane of this host under this assay?

Allowed labels: FACS bright/dim, in-cell GFP, in-gel GFP above a source-defined threshold, TargetTrack status `expressed` when the center used a membrane-protein protocol.

Not allowed as A: purification yield, SEC monodispersity, PDB deposition, “soluble fraction” of a WRAP fusion.

### B — produced as studyable material

Did a production pipeline report that the protein was extracted and purified (optionally crystallized)?

Primary public proxy: TargetTrack statuses `soluble` / `purified` / `crystallized` / `in PDB`, always conditioned on the target having been attempted (`cloned` or later). The estimand is **P(purified | expressed)** or **P(purified | cloned)**, reported separately, never as a single “success” bit mixed with A.

Crude FSEC / FSEC-TS is what a new 1000-plex campaign would measure. We do not invent FSEC traces. TargetTrack is not FSEC.

### B-conditions — how it was purified, given that it was

Chromatography buffer recipes from PurificationDB (when a dump exists) plus GPCRdb construct annotations (solubilization detergent, chromatography, crystallization chemicals) and four published detergent screens (Lantez 2015, Lin 2016, Högbom/Sjöstrand 2017, Kotov 2019). Every GPCRdb row is a solved-GPCR success. The screens are small-n and/or success-biased. **None of these tables train “will it purify.”** They are a starting-recipe prior after the playbook has already chosen `standard` vs WRAP vs redesign. Högbom Table 2 per-protein FSEC grades are figures, not a public matrix — we keep the 60-protein list, the 16-detergent list, and the family-level claims from the text.

### C — historically solved

Is there a membrane-embedded structure in PDBTM or mpstruc?

This is a popularity- and tractability-biased terminal label. It is a contrast with UniProt Swiss-Prot TM proteins, not a purification oracle.

---

## 2. Units

- **Construct** — one amino-acid sequence as assayed (fusion tags stripped when the source provides the insert sequence; otherwise recorded as `tag_unknown`).
- **Target** — one TargetTrack target ID; a target may have several trials.
- **Status** — the latest reported TargetTrack status for that target, mapped onto the ordered ladder in Section 3.
- **Architecture** — `bitopic_helical` | `polytopic_helical` | `beta_barrel` | `monotopic` | `unknown`, from UniProt/TOPDB/PDBTM when present, else from Kyte–Doolittle TM count.

---

## 3. TargetTrack status mapping

Ladder, low to high (a later status implies the earlier ones were attempted; it does not imply they succeeded in every trial):

1. `selected`
2. `cloned`
3. `expressed` → question A (center-dependent)
4. `soluble` → weak B
5. `purified` → B
6. `crystallized` / `diffraction` / `NMR` → B plus structural intent
7. `in PDB` / `structures` → C-overlapping, still not pooled with C from PDBTM

Membrane filter (either is sufficient):

- Center name matches a membrane-protein effort (NYCOMPS and any center whose deposited protocol text contains membrane-protein production), **or**
- Sequence has ≥2 predicted TM helices (Kyte–Doolittle, window 19, threshold 1.6), **or**
- UniProt accession joins to `ft_transmem`.

Proteins with 1 predicted TM are kept in a separate `bitopic` slice and not mixed into the polytopic B model.

---

## 4. Sources and what they may answer

| Source | Question | Host / bias |
|---|---|---|
| Curnow labelled FACS | A | E. coli, one designed scaffold, selected tails |
| Curnow unlabelled remainder | A inference only | same scaffold |
| TargetTrack TM-filtered | A (expressed), B (purified \| cloned) | mixed centers, mostly prokaryotic |
| Daley / Hammon / Drew GFP | A transfer | E. coli or yeast, small n |
| PurificationDB | B-conditions | human PDB, success-only, detergent NER incomplete |
| GPCRdb constructs | B-conditions | solved GPCRs, success-only solubilization/xtal recipes |
| Högbom 2017 / Kotov 2019 / Lin 2016 / Lantez 2015 | B-conditions | small published screens; not purify-vs-fail |
| PDBTM / mpstruc | C | solved structures |
| UniProt Swiss-Prot TM | C background; ranking universe | reviewed TM, all organisms |
| TOPDB / HTP | topology features | experimental + predicted topology |

A row always carries `source` and `question`. Duplicate sequences across sources are kept as separate rows.

---

## 5. Features (v1)

Computed only from amino-acid sequence plus topology annotation or Kyte–Doolittle:

- length, AA composition (20), mean Kyte–Doolittle, net charge (pH 7)
- predicted TM count, TM residue fraction, longest TM, longest extra-membrane loop
- N- and C-terminal segment orientation (GFP-fusion placement)
- N-X-S/T sequon count (not in predicted TM)
- Cys count

No embeddings, no DeepTMHMM-as-a-service, no codon features unless a nucleotide sequence is supplied by the source (Curnow DNA is optional and not used in the primary A model).

**Amendment (embeddings):** ESM-2 650M mean-pooled vectors may be concatenated to v1 features or used alone. They are reported under the same splits. A rise in random-split AUC without a rise in center/organism AUC is leakage, not a better HTS prior.

---

## 6. Splits

Random split is reported as an **optimistic** control. Primary claims use:

- **Neighborhood / Hamming cluster** on Curnow (equal-length designed library).
- **Center** on TargetTrack (NYCOMPS vs other centers).
- **Organism domain** (prokaryote vs eukaryote) on natural-protein tables.
- **Architecture** (helical polytopic vs beta barrel vs bitopic).
- **Time** for C (PDB deposition year), when the year is present.

Curnow-trained scores applied to Swiss-Prot TM or TargetTrack are labeled **out of distribution**. They are not “predicted yield.”

---

## 7. Models and metrics

- Logistic regression and random forest, CPU, scikit-learn defaults except `n_estimators=200` for RF and class-balanced weights.
- Metric: ROC-AUC. Also report PR-AUC and the positive rate in train/test.
- Sample size floor: a split is skipped if either side has <30 positives or <30 negatives.
- Question A, B, and C models are separate objects. No multi-task head in v1.

---

## 8. WRAP rescue (not a claim)

Baker WRAPs are genetically fused amphipathic proteins that solubilize a TM target in the E. coli cytoplasm without detergent. v1 scores **amenability**, it does not design WRAPs.

Heuristic inputs: predicted TM belt hydrophobicity, n_TM ≥ 1, length, whether at least one terminus is extra-membrane or cytoplasmic (fusion-able), architecture. Ranked list of B-failures and Swiss-Prot TM proteins. Vocabulary: `wrap_amenable` is allowed; `solubilized`, `designed WRAP`, and `detergent-free structure` are not.

---

## 9. Stop rule (Gate 0)

Before fitting:

1. Print n per source and question.
2. For B: n of TM-filtered TargetTrack targets that reached `cloned` and n that reached `purified`. If purified positives <100 **or** they come from a single center, **do not claim a B model.** Report the census and stop that arm.
3. For B-conditions: n of PurificationDB UniProts that join to Swiss-Prot TM or PDBTM, and n whose buffer lists a membrane detergent (DDM, DM, OG, LDAO, LMNG, GDN, FC12, CHAPS, tritons, cholate). If the TM slice is <50, keep PurificationDB as a lookup table, not a model. GPCRdb unique PDBs and literature screens are a **recipe prior** when n_PDB ≥ 50; they still do not authorize a purify-vs-fail classifier.
4. For A: Curnow labelled n is already known to be ~2e3 on one scaffold — sufficient for a local A model, insufficient for transfer. Transfer tests are reported as transfer, not as A skill.
5. For C: require ≥500 mpstruc/PDBTM unique proteins and ≥2,000 Swiss-Prot TM background proteins.

---

## 10. What a positive result is allowed to mean

- A model of Curnow FACS predicts FACS on held-out Curnow neighborhoods — local expression on that scaffold.
- A drop from random-split AUC to center/family/organism AUC is evidence of leakage, not of a worse model.
- P(purified | cloned) on TargetTrack TM is “this center’s pipeline recovered protein,” not “FSEC-monodisperse” and not “functionally folded.”
- C vs Swiss-Prot TM is “resembles proteins that have been solved,” not “will express in my host.”
- The playbook (B-action) recommends a next experiment class: `standard` vs `wrap_rescue` vs `redesign` vs `deprioritize`. It is not a buffer recipe and not a guarantee of purified protein.
- The detergent prior attached to a playbook row is a starting extract/SEC recipe (DDM ± CHS; LMNG-class rescue; avoid fos-choline/PEG as stability detergents). GPCRdb recipes are GPCRs. Högbom/Kotov/Lin/Lantez are small published screens, not a 1000-protein FSEC table.

---

## Amendments

- 2026-09-11 Gate 0: UniTmp PDBTM/TOPDB bulk XML and PurificationDB dumps were unreachable at ingest time. Question C uses mpstruc (4,285 unique PDB IDs) joined to Swiss-Prot TM. B-conditions stays lookup-only (n=0). GFP-paper supplements were not present as machine-readable tables.
- 2026-09-11 B-action playbook: `mpx-playbook` trains P(purified | expressed) on TargetTrack membrane centers and layers WRAP amenability plus construct levers. Allowed outputs are the four action tokens above. Forbidden: claiming detergent conditions, designed WRAPs, or FSEC success.
- 2026-09-11 embeddings: ESM-2 650M mean-pool on the unique-sequence universe (Swiss-Prot TM + Curnow + TOPDB + TargetTrack). Same leakage splits as composition features. Modes: `compose`, `esm`, `both`. Embeddings are not a new question and are not allowed to be reported as A skill on natural proteins unless the GFP transfer split is shown.
- 2026-09-11 GFP transfer tables: Daley 2005 Table S1 GFP fusions with GFP/ml > 0 (n=579, median split); Hammon 2009 Supporting Table 1 (313/314 parsed; label FSU ≥ 60,000 as in the paper, 64 positives). Sequences joined to Swiss-Prot TM by E. coli gene (Daley) or fetched from UniProt/NCBI by accession (Hammon).
- 2026-09-11 TOPDB XML ingested from a local dump (unitmp bulk URLs still 404). PDBTM bulk XML remains missing.
- 2026-09-12 Hammon accessions: UniProt REST + NCBI protein FASTA filled **313/313** sequences (263 UniProt, 50 RefSeq). Daley remains a Swiss-Prot TM gene join (**428/579** with sequence). These rows are question A **transfer** only.
- 2026-09-12 embeddings: ESM-2 650M mean-pool universe n=101807 (Swiss-Prot TM + Curnow + TOPDB + TargetTrack), truncated 1022 aa. Modes `compose` / `esm` / `both` on the same leakage splits. Recorded RF holdout: A cluster compose **0.868** vs ESM **0.835**; B center compose **0.596** vs ESM **0.670** vs both **0.681**; C organism compose **0.552** vs ESM **0.603**. Curnow → Daley best **0.570**; Curnow → Hammon compose RF **0.441** (ESM on 300/313 sequences). Transfer is not A skill.
