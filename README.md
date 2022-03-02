# Run on farm

1. Create your yaml following [this](https://gitlab.internal.sanger.ac.uk/tl10/workflow-registration/-/blob/lsf/template.yaml)
2. Log into an interactive session `bsub -q imaging -n1 -M4000 -R"select[mem>4000] rusage[mem=4000]" -Is $SHELL`
2. Run with `/lustre/scratch117/cellgen/team283/tl10/workflow-registration/run.sh [YOUR_YAML]`

