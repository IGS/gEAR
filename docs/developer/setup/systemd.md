# systemd

gEAR's background services run as systemd units. The unit files are in `systemd/` in the repository and are copied to `/etc/systemd/system/` on each server.

This page covers the general systemd ideas. The steps for each gEAR service are in that service's doc:

| Units | Doc |
| --- | --- |
| RabbitMQ consumers: `gear-consumers.target`, `*-consumer.target`, `*-consumer@.service`, `gear-consumers.slice` | [RabbitMQ Consumers](../services/rabbitmq_consumers.md): install, start/stop, [rolling out updates](../services/rabbitmq_consumers.md#rolling-out-updates), [changing unit files](../services/rabbitmq_consumers.md#changing-unit-files), scaling |
| Spatial Panel: `spatial-panel.service` | [Spatial Panel](../services/spatial.md) |
| Consumer memory limits (`gear-consumers.slice`, per-service drop-ins) | [RabbitMQ setup](./rabbitmq.md#preventing-runaway-memory) |
| Apache override | [Apache](./apache.md) |

## Templates and targets

- A unit whose file name has an `@` (for example `worker@.service`) is a **template**. Each instance gets a name with a number after the `@`, such as `worker@1.service` or `worker@2.service`, and runs as its own process. `%i` in the file is replaced with that number. See <https://www.stevenrombauts.be/2019/01/run-multiple-instances-of-the-same-systemd-unit/>.
- A **target** groups other units so they can be handled with one command.

A target and its members are linked in two directions, and each link controls different commands:

- `Wants=` on the **target** lists its members. It is used only when starting: starting the target starts any member that is not running. It does nothing for stop or restart.
- `PartOf=` on the **member** names the target. It is used only for stop and restart: stopping or restarting the target also stops or restarts the member.

A target that should control its members fully needs both. The gEAR consumer units have both at every level.

## Common commands

| Command | What it does |
| --- | --- |
| `sudo systemctl start <unit>` | Start now. Does nothing to a unit that is already running, so it does **not** load new code. |
| `sudo systemctl restart <unit>` | Stop, then start. Use this to load new code. |
| `sudo systemctl stop <unit>` | Stop now. |
| `sudo systemctl enable <unit>` | Start at boot. Does not start it now; also run `start`, or use `enable --now`. |
| `sudo systemctl daemon-reload` | Read the unit files again after copying or editing them. **Running units keep their old settings until they are restarted.** |
| `sudo systemctl reset-failed <unit>` | Clear a `failed` state (for example `start-limit-hit`) so the unit can be started again. |
| `systemctl status <unit>` | Show state, start time, and recent log lines. Accepts patterns such as `'projectr-consumer@*'`. |
| `sudo journalctl -u <unit>` | Full logs. Add `-f` to follow, `-n 100` for the last 100 lines. |

## Notes

- For a long-running process such as a RabbitMQ consumer, use `Type=simple` in `[Service]`. With `Type=forking`, systemd waits for the process to fork and exit, so `systemctl start` never returns.
- When changing a unit file, stop the running unit **before** renaming or deleting its file. Otherwise systemd shows the unit as `not-found` while the process keeps running.
