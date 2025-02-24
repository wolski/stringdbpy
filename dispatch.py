# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "polars",
#     "bfabric-app-runner",
# ]
#
# [tool.uv.sources]
# bfabric-app-runner = { git = "https://github.com/fgcz/bfabricPy.git", subdirectory = "bfabric_app_runner", rev = "main" }
# ///
import polars as pl
from bfabric_app_runner.dispatch.dispatch_resource_flow import ResourceDispatcherCLI


def dispatch_strategy(resources_df, workunit_definition):
    resources_df = resources_df.filter(pl.col("filename").str.ends_with(".zip"))
    extra_inputs = [
        {
            "type": "static_yaml",
            "filename": "params.yml",
            "data": workunit_definition.execution.raw_parameters,
        }
    ]
    return resources_df, extra_inputs


if __name__ == "__main__":
    ResourceDispatcherCLI(dispatch_strategy).run()
