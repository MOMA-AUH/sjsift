# README screenshot

This example creates the [README screenshot](../../docs/images/splice-view.png) from a real offline HTML report.
It uses METex14, the exon 13-to-15 junction, and both reference junctions from the repository's GRCh38 catalog.
The read names, sequences, alignments, and STAR-format support counts are synthetic.
The counts are invented examples, not results from a STAR run.

The report contains six variant alignment records, 60 same-donor reference records, and 240 same-acceptor reference records.
One shorter same-donor record is multimapping. All other records are unique.
Reference reads have varied starts, ends, and lengths.
The synthetic STAR-format counts match these totals, so the larger reference counts produce visibly thicker junction arcs.
The initial view shows all six variant reads and ten reads per reference group. The remaining reads are available through pagination.
The example does not represent a biological sample or clinical validation.

## Recreate the image

Use the Python development environment and Node.js 24 described in the [development guide](../../docs/DEVELOPMENT.md#local-setup).
From the repository root, install the existing browser tools:

```bash
npm ci --prefix tests
tests/node_modules/.bin/playwright install --with-deps chromium
```

Generate the report in a new directory:

```bash
python examples/readme-preview/generate.py /tmp/sjsift-readme-preview
```

Open `/tmp/sjsift-readme-preview/synthetic.html` in a browser and select `METex14` to inspect the interactive report.
The directory also contains the synthetic inputs and output tables.
The generator refuses to reuse an existing directory.

Capture the variant detail panel:

```bash
node examples/readme-preview/capture.cjs \
  /tmp/sjsift-readme-preview/synthetic.html docs/images/splice-view.png
```

The capture command replaces the PNG at the specified path.
It uses a 1200-pixel browser viewport and twice the pixel density for readable text.
It preserves the report's content and styles, checks the visible reads and arcs, and rejects browser errors or external requests.
Chromium and system fonts can affect the rendered pixels.
Keep the generated HTML, BAM, index, and TSV files outside the repository.
