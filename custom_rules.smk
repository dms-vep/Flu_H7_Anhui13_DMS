"""Custom snakemake rules for additional analyses."""

# Read the wrapped heatmap config
wrapped_heatmap_config_file = "data/wrapped_heatmap_config.yaml"
if os.path.isfile(wrapped_heatmap_config_file):
    with open(wrapped_heatmap_config_file) as f:
        wrapped_heatmap_config = {
            key: val
            for (key, val) in yaml.YAML(typ="safe", pure=True).load(f).items()
        }
else:
    wrapped_heatmap_config = {}


rule wrapped_heatmap:
    """Create wrapped heatmap visualization."""
    input:
        data_csv=lambda wc: wrapped_heatmap_config[wc.wrapped_hm]["data_csv"],
    output:
        chart_html="results/wrapped_heatmaps/{wrapped_hm}.html",
    params:
        params_dict=lambda wc: wrapped_heatmap_config[wc.wrapped_hm],
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    log:
        notebook="results/notebooks/wrapped_heatmap_{wrapped_hm}.ipynb",
    notebook:
        "notebooks/wrapped_heatmap.py.ipynb"


# Add wrapped heatmaps to docs
if wrapped_heatmap_config:
    docs["Wrapped heatmaps"] = {
        "Heatmap HTMLs": {
            wrapped_hm: rules.wrapped_heatmap.output.chart_html.format(wrapped_hm=wrapped_hm)
            for wrapped_hm in wrapped_heatmap_config
        },
    }
