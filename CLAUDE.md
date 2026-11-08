# MitoM - notes for agents

- **Tests: check `knowledge/TEST_REGISTRY.md` before running any benchmark / real-data test** (two machines run tests:
  `deyuan` = /media/deyuan/18T_B/yzz/mitom_work, `cody` = /media/cody/18T_R1). Do not repeat a test that is listed
  unless the code it measured has changed; say which row you rely on. After any test add one row there (date, machine,
  question, data, tools/versions, key numbers, conclusion/commit, where recorded) and keep details in
  `benchmark/results/<name_date>/README.md`. Commit both so the other machine sees them.
- Run errors: cause and fix go to the troubleshooting tables (docs/USAGE_zh.md section 9, docs/MITOS2.md section 6).
- Docs for users: Chinese + English, short and plain. A measured number (accuracy, time, memory, counts from a benchmark) is written once, in a table with bilingual headers. The Chinese paragraph and the English paragraph point at that table and do not repeat the number. Option defaults and rule thresholds stay next to the option or the rule, in one place. Full options: `python3 scripts/make_options_doc.py` after CLI changes
  (tests/test_docs.py checks it).
