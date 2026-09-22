# Use TOML for variant-definition catalogs

Variant-definition catalogs use TOML because they combine catalog-level metadata with repeated, typed variant entries and must remain straightforward for humans to review and edit. TOML avoids a custom TSV metadata convention and YAML's additional parser and implicit typing, while Python 3.11 and newer can parse it from the standard library.
