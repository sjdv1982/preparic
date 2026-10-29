"""Subprocess entry point for the W1 service-free local worker acceptance test."""

import asyncio
import json
import os


def main():
    import seamless
    import seamless.config
    import seamless_config
    from seamless.transformer import delayed
    from seamless_config.extern_clients import collect_remote_clients
    from seamless_transformer import worker

    def forbidden_init(*args, **kwargs):
        raise AssertionError("local execution must not initialize configured services")

    seamless.config.init = forbidden_init
    seamless_config.init = forbidden_init

    import backend

    @delayed
    def child_identity(index):
        import os
        import socket
        import time

        time.sleep(0.1)
        return {
            "index": index,
            "pid": os.getpid(),
            "tcp_disabled": socket.socket.connect.__module__ == "sitecustomize",
        }

    async def run_tasks():
        return await asyncio.gather(*(child_identity(i).task() for i in range(16)))

    try:
        configured = backend.configure_local_execution()
        handles = worker._require_manager()._handles
        assert configured == 8
        assert len(handles) == 8
        assert backend.configure_local_execution() == 8
        assert len(worker._require_manager()._handles) == 8
        results = asyncio.run(run_tasks())
        clients = collect_remote_clients("w1-local-spawn")
        assert clients == {"database": [], "buffer": []}, clients
        print("W1_SPAWN_RESULT=" + json.dumps({
            "configured_workers": configured,
            "worker_handles": len(handles),
            "parent_pid": os.getpid(),
            "results": results,
            "clients": clients,
            "persistent_cache": seamless_config.get_seamless_cache(),
        }), flush=True)
    finally:
        seamless.close()


if __name__ == "__main__":
    main()
