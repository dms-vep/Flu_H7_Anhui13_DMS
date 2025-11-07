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
        data_csv=lambda wc: wrapped_heatmap_config[wc.heatmap]["data_csv"],
        nb="scratch_notebooks/wrapped_heatmap.py.ipynb",
    output:
        chart_html="results/wrapped_heatmaps/{heatmap}.html",
        nb="results/notebooks/wrapped_heatmap_{heatmap}.ipynb",
    params:
        params_dict=lambda wc: wrapped_heatmap_config[wc.heatmap],
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    log:
        "results/logs/wrapped_heatmap_{heatmap}.txt",
    run:
        import yaml as run_yaml
        # Create parameters for papermill
        pm_params = {
            "data_csv": str(input.data_csv),
            "chart_html": str(output.chart_html),
        }
        # Add all parameters from config
        pm_params.update(params.params_dict)

        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            run_yaml.dump(pm_params, f)
            yaml_file = f.name

        shell(
            f"papermill {{input.nb}} {{output.nb}} "
            f"-f {yaml_file} "
            f"&> {{log}} && rm {yaml_file}"
        )


# Add wrapped heatmaps to docs
if wrapped_heatmap_config:
    docs["Wrapped heatmaps"] = {
        heatmap.replace("_", " ").replace("/", " - "): rules.wrapped_heatmap.output.chart_html.format(
            heatmap=heatmap
        )
        for heatmap in wrapped_heatmap_config
    }

    # Add to other_target_files so they get built
    other_target_files.extend([
        rules.wrapped_heatmap.output.chart_html.format(heatmap=heatmap)
        for heatmap in wrapped_heatmap_config
    ])
