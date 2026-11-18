# Current membership and split rules

## S1 Data construction and source reconciliation

We organized ProNAB affinity records by molecular construct, source study and experimental context. Absolute cohort A comprised reported dissociation constants, Kd, with usable units and recorded construct identity. Source-specific measurement or estimation procedures are retained in the landscape correspondence table; for example, the EF-Tu records describe Kd inferred from dissociation rates under a fixed association-rate assumption. Affinity was expressed as pKd = −log10(Kd/M). Repeated measurements assigned to the same unit were summarized by their median on the pKd scale; distinct constructs and conditions retained separate identities.

Protein-mutation cohort B linked each single substitution to a WT anchor from the same study and experimental conditions. Multi-point substitutions, simultaneous changes on both molecular sides and unresolved WT links were excluded. RNA mutants could enter absolute evaluation. A landscape denotes one linked WT–mutant experimental system. B-shape required at least five distinct mutated protein positions both in the eligible landscape and after restricting to shared prediction coverage. Smaller landscapes remained in ordinary mutation-effect evaluation.

Structural inputs used the mapped protein and RNA chains. In the controlled CoPRA comparison, sequence substitutions retained WT coordinates; the manifest placed each change on the first compatible matched chain. The official parser selected the first coordinate model, and the input transform selected a 256-node interface patch. Other structural applications used their recorded mutant-building protocols. Each comparison retained predictions from all participating methods and, for OOF models, all three seeds. Linked absolute analyses counted each WT unit once.

Table S30 links the stable L01–L23 archive identifiers to original-source identifiers, experimental constructs and conditions, PDB chains and numbering maps. Its companion evidence index records parsed chain lengths, sequence differences and known input discrepancies. Reference sequence length is distinguished from assay fragment length, and multiple structure copies are distinguished from the experimentally mutated assembly. The principal CoPRA scope retains the 1JBR A-copy mutations; its source-reconstructed crop retains all eleven changed positions. Experimental mutated-copy stoichiometry remains unresolved, so the existing whole-landscape exclusion is reported alongside the retained scope (Table S31). Later structural applications treat their 1JBR inputs as separate position controls.

Source inspection identified seven censored Kd units in the frozen evaluation: six protein mutants and one RNA mutant whose lower bounds had been transcribed as point values (Dominguez et al. 2010). A subsequent source-figure check resolved L19/L20 as duplicate Hfq records of L11/L10 and retained the source-consistent WT anchors, removing ten duplicate pairs and twelve endpoint units (S1). The principal summaries therefore use 752 absolute units, 338 pairs and 220 B-shape pairs. The original fitted models and training manifests are unchanged; the corrections update evaluation membership and weights. The same Hfq alias rule is propagated within each affected method-specific comparison (Tables S3–S7).

**Table S1. Current principal evaluation populations.** All rows exclude source-confirmed censored observations and Hfq aliases. Linked absolute units deduplicate WT anchors.

| Population | Units or pairs | Landscapes | Families |
|---|---:|---:|---:|
| Absolute affinity | 752 units | — | 113 |
| Ordinary mutation effects | 338 pairs | 70 | 52 |
| Local mutation effects | 220 pairs | 21 | 20 |
| Absolute endpoints linked to ordinary pairs | 408 units | — | — |
| Absolute endpoints linked to local pairs | 241 units | — | — |

L10/L20 and L11/L19 each duplicate one set of seven-nucleotide OxyS-fragment FP measurements (Wang et al. 2015). The ten paired mutant values match the original figures and their corresponding alias values. The canonical L10/L11 WT values also match those figures; the L20/L19 WT anchors do not. All twelve corresponding CoPRA/additive endpoint predictions and consumed input identities are identical across the aliases. We retain source-consistent records and collapse aliases without averaging conflicting WT values.

The source describes full-length Hfq in FP and Hfq65 (residues 1–65) in crystallography. Its pH-7 phosphate/NaCl/EDTA/MgCl2 buffer belongs to a separate full-length OxyS_U8 EMSA experiment (Figures 6E/F). The seven-base FP measurements are in Figures 1B/C and 6B/C/D. The article refers to prior FP protocols rather than restating exact temperature and buffer. Different database method labels and condition fields therefore do not establish separate experiments here.

**Table S2. Hfq source identity and evaluation weight.** Each row has five substitutions; the two RNA fragments remain separate experimental systems.

| Archive IDs | Source experiment | Source figures | Retained record / weight |
|---|---|---|---|
| L10; L20 | A_us = AUAACUA, FP | 1B; 6B | L10 retained; L20 alias, zero additional weight |
| L11; L19 | A_ds = AACUAAA, FP | 1C; 6C/D | L11 retained; L19 alias, zero additional weight |

The source rule is propagated only within each method's available members. The two source-unique Hfq landscapes share their family's weight equally. Constant WT offsets disappear under centering, so family-macro centered error, ordering and response SD are unchanged; ordinary and endpoint metrics use the retained source-consistent WT anchors. PRA, M2, SVSBI, GBNL, ProRB and PANTHER contain no removed aliases in their current comparison scopes. Hfq occurs only in the training partition of the three-arm study, whose training-fit panels retain as-trained weights; its four-family local validation scope is unaffected.

The censored records originate from hnRNP F qRRM measurements (Dominguez et al. 2010). Lower-bound Kd values are retained as bounds in the source record and excluded from point-error evaluation. Family-fold fits retain their original training labels; source corrections alter evaluation membership and weights. The source index records both training history and final evaluation membership.


The import adapter executes the saved membership and exclusions; it does not infer experimental links from sequence similarity. Training membership remains the as-fitted 771-unit manifest. Evaluation changes do not retrain those fits. Exact row-level maps and fold assignments must accompany authorized input delivery.
