# Results

The `csv/` folder contains the final candidate and submitted CSV trajectories kept
from the project root.

## Best Preserved Candidates

The top candidates were renamed to include a short story tag while keeping the score:

| File | Score in filename | Story tag |
| --- | ---: | --- |
| `topscore_J251.569127_...csv` | 251.569127 | Best overall score from the final morning run |
| `runnerup_J249.609577_...csv` | 249.609577 | Second-best run with a near-identical backbone sequence |
| `third_J241.274564_...csv` | 241.274564 | Third-best score with a slightly different mid-tour ordering |
| `grandtour_J183.936289_...csv` | 183.936289 | Grand Tour sequence visiting all bodies (incl. Nyxar) |

Other preserved candidates keep their original naming (e.g. `best_J175.982794_...csv`
and `rank12_J233.055070_...csv`). Additional files named by event score, such as
`222p075764_events_26.csv`, are also kept because they were part of the final
submission set.

## Validator

The final score is assigned by the online validator after uploading a CSV file. The
repository contains the CSV generation and trajectory search code, but it does not
include the validator itself.
