"""Transform the MedQIP ABP narrative extract from long to wide format.

Input : one row per (assessment_event_id, question), with event-level
        metadata repeated on every row.
Output: one row per assessment_event_id, with each question's response
        spread into its own column.

For each question the response `value` becomes `<question>_value`; where the
question is answered on an anchored scale, the anchor text is also kept as
`<question>_anchor`. Instrument title/id can differ across questions within
one event (the narrative questions belong to the "General Feedback"
instrument), so the unique values are joined with " | ".

Usage: python long_to_wide.py <input.xls> <output_basename>
"""

import re
import sys

import pandas as pd

EVENT_ID = "assessment_event_id"

# Columns that are constant within an assessment event (verified in the data:
# each has exactly 1 unique value per event).
EVENT_COLS = [
    "program_id_anon",
    "activity_title",
    "activity_id",
    "activity_event_id",
    "creator_id_anon",
    "activity_event_date",
    "activity_event_created_date",
    "evidence_recorded_date",
    "assessment_event_creation_date",
    "expiration_date",
    "rater_id_anon",
    "rater_role_during_event",
    "subject_id_anon",
    "subject_role_during_event",
    "subject_pgy_now",
    "subject_pgy_during_event",
    "created_by",
]

# Question-level columns that vary within an event but can hold 2 values
# (main instrument + the narrative-only "General Feedback" instrument).
MULTI_VALUE_COLS = ["instrument_title", "instrument_id"]


def slug(title: str) -> str:
    """'Duration of Observation' -> 'duration_of_observation'."""
    return re.sub(r"[^a-z0-9]+", "_", title.lower()).strip("_")


def main(infile: str, outbase: str) -> None:
    df = pd.read_excel(infile)

    # Sanity check: the pivot key must be unique per event.
    dup = df.duplicated(subset=[EVENT_ID, "question_title"]).sum()
    if dup:
        raise ValueError(f"{dup} duplicate (event, question_title) rows; pivot would be ambiguous")

    event = df.groupby(EVENT_ID, sort=False)[EVENT_COLS].first()
    for col in MULTI_VALUE_COLS:
        event[col] = df.groupby(EVENT_ID, sort=False)[col].agg(
            lambda s: " | ".join(pd.unique(s.dropna().astype(str)))
        )

    values = df.pivot(index=EVENT_ID, columns="question_title", values="value")
    anchors = df.pivot(index=EVENT_ID, columns="question_title", values="answer_option_anchor")

    # Keep an anchor column only for questions that ever have anchor text
    # that differs from the value itself.
    keep_anchor = [
        q for q in anchors.columns
        if anchors[q].notna().any() and not anchors[q].dropna().astype(str).equals(
            values[q].reindex(anchors[q].dropna().index).astype(str)
        )
    ]
    anchors = anchors[keep_anchor]

    values.columns = [f"{slug(q)}_value" for q in values.columns]
    anchors.columns = [f"{slug(q)}_anchor" for q in anchors.columns]

    # Interleave value/anchor columns so each question's pair sits together.
    ordered = []
    for vcol in values.columns:
        ordered.append(values[vcol])
        acol = vcol[: -len("_value")] + "_anchor"
        if acol in anchors.columns:
            ordered.append(anchors[acol])

    wide = pd.concat([event] + ordered, axis=1).reset_index()
    # reset_index puts the id first; move it back after the event metadata? No —
    # keep assessment_event_id as the leading identifier column.

    wide.to_csv(f"{outbase}.csv", index=False)
    wide.to_excel(f"{outbase}.xlsx", index=False)
    print(f"Long : {df.shape[0]} rows x {df.shape[1]} cols")
    print(f"Wide : {wide.shape[0]} rows x {wide.shape[1]} cols")
    print(f"Wrote {outbase}.csv and {outbase}.xlsx")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
