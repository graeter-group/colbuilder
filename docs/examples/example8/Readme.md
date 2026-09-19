<div align="center">
    <h1>Example 8</h1>
    <p>Bringing your own sequence: converting a PDB to FASTA with <code>pdb2fasta</code> for a custom species.</p>
</div>

Demonstrates the "custom species" path: providing your own `fasta_file` for a
species name colbuilder doesn't recognize, using the bundled `pdb2fasta`
console script to produce that FASTA from an existing PDB.

## Usage

```bash
# 1) Build a plain (non-crosslinked) reference structure to convert
colbuilder --config_file config_step1_base_sequence.yaml

# 2) Convert it to FASTA with pdb2fasta (a colbuilder console script)
pdb2fasta rattusnorvegicus_N_NONE_C_NONE.pdb > custom_species.fasta

# 3) Feed that FASTA back in under a species name colbuilder doesn't know,
#    to build a fresh structure for it
colbuilder --config_file config_step2_custom_species_sequence.yaml
```

## Files

- `config_step1_base_sequence.yaml` — `crosslink` unset, sequence generation for
  `rattus_norvegicus`, just to produce a plain PDB to convert (any PDB with standard
  CA-atom records works here — this isn't rattus_norvegicus-specific).
- `config_step2_custom_species_sequence.yaml` — `species: "my_custom_species"` (not in colbuilder's
  built-in species list) with `fasta_file: "custom_species.fasta"` set explicitly. Without
  `fasta_file`, this fails at config validation ("Must provide fasta_file when using
  custom species").

## Output Files

- `rattusnorvegicus_N_NONE_C_NONE.pdb` — the plain source structure.
- `custom_species.fasta` — produced by `pdb2fasta`; three chain entries
  (`>rattusnorvegicus_N_NONE_C_NONE.pdb:A/B/C`), standard residues mapped to one-letter
  codes (hydroxyproline `HYP` → `O`, per `pdb2fasta`'s table).
- `custom_species_alignment.fasta`, `custom_species_N_NONE_C_NONE.pdb` — sequence-generation
  output for the custom species, named from the `species` value in step 2's config.

## Notes

- This example doesn't apply crosslinks to the custom species: `n_term_type`/`c_term_type`
  combinations are looked up by exact species-name match in `crosslinks.csv`, which raises
  `"No crosslinks found for species: ..."` for any name not already a row there. Crosslinking
  a genuinely novel species means also adding matching rows yourself, with residue positions
  specific to your own sequence — a separate, more advanced task this example doesn't cover.
- `pdb2fasta`'s residue-mapping table only handles standard amino acids plus colbuilder's
  known crosslink marker codes; give it a PDB with already-applied, unrecognized modified
  residues and those positions map to `-` (gap) in the output FASTA.
