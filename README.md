# Nanome - Data Table

A Nanome plugin to view multi-frame structure metadata in a table

## Dependencies

[Docker](https://docs.docker.com/get-docker/)

## Usage

To run Data Table in a Docker container:

```sh
$ cd docker
$ ./build.sh
$ ./deploy.sh [run_args]
```

### Optional arguments

- `--https`

  Enable HTTPS using a self-signed certificate. If port is not set, port will default to 443.

- `-u url` or `--url url`

  The url opened when the plugin is run to access the Web UI. Example: `-u table.example.com`

- `-w port` or `--web-port port`

  The port to use for the Web UI. Example: `-w 8080`

  Some OSes prevent the default port `80` from being used without elevated permissions, so this option may be used to change to an allowed port.

## Development
To run Data Table with autoreload:

```sh
$ python3 -m pip install -r requirements.txt
$ python3 run.py -r [run_args]
```

## License

MIT
