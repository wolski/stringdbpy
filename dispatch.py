# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "bfabric",
#     "bfabric-app-runner",
#     "cyclopts",
#     "pyyaml",
# ]
# ///
from bfabric import Bfabric
from bfabric.entities import Resource
from bfabric.experimental.workunit_definition import WorkunitDefinition
from bfabric_app_runner.dispatch.dispatch_single_resource_flow import (
    DispatchSingleResourceFlow,
    ConfigDispatchSingleResourceFlow,
)
from pathlib import Path
import yaml
import cyclopts


class Dispatch(DispatchSingleResourceFlow):
    def dispatch_job(
        self, resources: list[Resource], workunit: WorkunitDefinition
    ) -> Path:
        chunk_dir = self._out_dir / "work"
        chunk_dir.mkdir(exist_ok=True, parents=True)
        inputs = [
            {
                "type": "bfabric_resource",
                "id": resource["id"],
                "filename": Path(resource["relativepath"]).name,
            }
            for resource in resources
        ]
        with (chunk_dir / "inputs.yml").open("w") as file:
            yaml.safe_dump({"inputs": inputs}, file)
        return chunk_dir


app = cyclopts.App()


@app.default
def main(workunit_ref: Path | int, out_dir: Path):
    client = Bfabric.from_config()
    config = ConfigDispatchSingleResourceFlow(filter_suffix=".zip")
    workunit_definition = WorkunitDefinition.from_ref(workunit_ref, client=client)
    dispatch = Dispatch(client, config, out_dir)
    dispatch.dispatch_workunit(workunit_definition)


if __name__ == "__main__":
    app()
