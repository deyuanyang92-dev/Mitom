# MitoM - notes for agents

- **Tests: check `knowledge/TEST_REGISTRY.md` before running any benchmark / real-data test** (two machines run tests:
  `deyuan` = /media/deyuan/18T_B/yzz/mitom_work, `cody` = /media/cody/18T_R1). Do not repeat a test that is listed
  unless the code it measured has changed; say which row you rely on. After any test add one row there (date, machine,
  question, data, tools/versions, key numbers, conclusion/commit, where recorded) and keep details in
  `benchmark/results/<name_date>/README.md`. Commit both so the other machine sees them.
- Run errors: cause and fix go to the troubleshooting tables (docs/USAGE_zh.md section 9, docs/MITOS2.md section 6).
- Docs for users: Chinese + English, short and plain. Full options: `python3 scripts/make_options_doc.py` after CLI changes
  (tests/test_docs.py checks it).
