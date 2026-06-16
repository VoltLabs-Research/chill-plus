# Chill+ Ice & Hydrate Classifier

Identifies hexagonal ice, cubic ice, interfacial ice, hydrate and interfacial hydrate via Chill+ bond-order correlation (Nguyen & Molinero 2015).

## Install

```bash
vpm install @voltlabs/chill-plus
```

## CLI

```bash
chill-plus <input_dump> [output_base] [options]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `<input_dump>` | yes | — | Input LAMMPS dump. |
| `[output_base]` | no | derived from input | Base path for output files. |
| `--cutoff` | no | `3.5` | O-O neighbor cutoff in angstroms. |

## Exports

| Output file | Exposure | Exporter → artifact |
|---|---|---|
| `{output_base}_atoms.parquet` | Chill+ Structure | AtomisticExporter → glb |
| `{output_base}_chill_plus.parquet` | Chill+ Summary | — (listing-only) |

---

Full input contract and examples: https://docs.voltcloud.dev/docs/plugins
