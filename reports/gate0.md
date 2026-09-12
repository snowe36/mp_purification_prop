# Gate 0 census

| Source | n | with sequence |
|---|---:|---:|
| curnow | 2176 | 2176 |
| curnow_unlabelled | 11503 | 11503 |
| gfp | 892 | 741 |
| uniprot | 80943 | 80943 |
| mpstruc | 4285 | 0 |
| pdbtm | 0 | 0 |
| topdb | 9558 | 9558 |
| targettrack | 20583 | 20583 |
| purificationdb | 0 | 0 |
| gpcrdb | 1289 | 0 |
| screens | 206 | 0 |

## B (TargetTrack)
- cloned 6747, expressed 5308, purified 2299
- funnel: {'selected': 20363, 'cloned': 6747, 'expressed': 5308, 'soluble': 3364, 'purified': 2299, 'in pdb': 83}
- claim B: **True**

## B-conditions (PurificationDB)
- n 0, TM slice 0, detergent rows 0

## B-conditions (GPCRdb + screens)
- GPCRdb rows 1289, solubilization 442, unique PDBs 268
- top detergents: {'DDM': 155, 'DM': 69, 'LMNG': 20, 'OG': 11, 'NG': 8, 'HTG': 4, 'sodium cholate': 3, 'digitonin': 1}
- top additives: {'CHS': 150, 'CHAPS': 6, 'sodium cholate': 4, 'Zn2+': 3}
- extract combos: {'DDM+CHS': 122, 'DM': 63, 'DDM': 29, 'LMNG+CHS': 19}
- screens: {'hogbom2017': {'targets': 60, 'detergents': 16, 'findings': 7}, 'kotov2019': {'targets': 9, 'detergents': 94, 'findings': 6}, 'lantez2015': {'targets': 0, 'detergents': 0, 'findings': 6}, 'lin2016': {'targets': 0, 'detergents': 0, 'findings': 8}}
- recipe prior (not a purify classifier): **True**

## C
- PDB ids 4285, Swiss-Prot TM 80943, join positives 2200

## GFP transfer tables (A, not pooled with Curnow)
- Daley n 579, joined 428
- Hammon n 313, joined 313

## TOPDB
- n 9558, with sequence 9558, PDB xref 7579

## Stop
- B: fit
- B-conditions: prior
- C: fit
- A: fit_local
