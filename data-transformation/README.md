# MedQIP ABP narrative extract: long → wide

Reshapes `medqip_research_data_extract_ABP_narrative_data_3.2026_deidentified.xls`
from long format (one row per question response, 2,264 rows × 28 columns) to wide
format (one row per `assessment_event_id`, 501 rows × 50 columns).

## Files

- `long_to_wide.py` — transformation script
- `medqip_ABP_narrative_data_wide.csv` / `.xlsx` — wide-format output

## How the reshape works

- **Row key:** `assessment_event_id` (501 unique events, 3–6 question rows each).
- **Event-level columns** (constant within an event — program, activity, dates,
  rater, subject, PGY, etc.) are kept once per row.
- **`instrument_title` / `instrument_id`** can hold two values per event (the
  main instrument plus "General Feedback (narrative questions only)"), so the
  unique values are joined with `" | "`.
- **Question columns:** each `question_title` becomes `<question>_value`
  (e.g. `feedback_value`, `supervision_value`). For questions answered on an
  anchored scale, the anchor text is kept alongside as `<question>_anchor`
  (e.g. `supervision_anchor`). `question_title` is unique within every event,
  so the pivot is unambiguous.

## Reproduce

```sh
python long_to_wide.py <input.xls> medqip_ABP_narrative_data_wide
```

Requires `pandas`, `xlrd` (for .xls input), and `openpyxl` (for .xlsx output).
