1. **HIGH — Import cleanup can delete the original input directory.** [mitom_pipeline.py:370](/media/cody/18T_R1/yzz/Mitom/mitom_pipeline.py:370)  
   `prepare_import()` replaces `args.input` with `00.imported/imported.gb`, but never saves `input_given`. Consequently, the new cleanup protects only the converted copy.  
   **Failing scenario:** Run `mitom import -i out/original_annotation -o out` against an existing MitoM output. Cleanup recursively deletes the original input directory, including auxiliary files that were never imported. Confirmed with mocked filesystem operations.  
   **Fix:** Save the original input before `prepare_import()` changes it, and carry that protection through all cleanup stages.

2. **MED — Linear tRNA filtering still accepts predictions spanning the circular origin.** [mito_trna.py:822](/media/cody/18T_R1/yzz/Mitom/mito_trna.py:822)  
   The linear branch checks only `start <= end`. MitoM-ML already represents origin-spanning predictions with `end > sequence_length`, so these pass the filter.  
   **Failing scenario:** For a 16,000-base FASTA labelled `topology=linear`, a MitoM-ML prediction at `15970..16035` survives and becomes `join(15970..16000,1..35)`, incorrectly joining the physical ends. Reproduced through the prediction adapter and feature builder.  
   **Fix:** Require `1 <= start <= end <= L` for linear records, and pass linear topology to predictors that support it.

3. **MED — Validation timeouts are reported as a successful submission step.** [mitom_steps.py:165](/media/cody/18T_R1/yzz/Mitom/mitom_steps.py:165)  
   The new `unchecked` counter recognizes only `MITOM.Table2asnFailed`; `validate_record()` returns `MITOM.Table2asnTimeout` on timeout.  
   **Failing scenario:** If table2asn times out on every record, `run_submission()` still returns `ok=True` and describes the records as checked. The error table contains timeout messages, but the step status incorrectly indicates success. Confirmed with an injected timeout result.  
   **Fix:** Count both failure codes as unvalidated, preferably using an explicit validation-completed status.

**Suggestions**

- Test input preservation through the complete import workflow, not just `clear_previous_work()`.
- Add coordinate tests covering both wrapped representations: `start > end` and `end > L`.