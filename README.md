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

## External Properties

This plugin supports configuration to fetch chemical properties from external sources. To get started, rename the `config.example.json` to `config.json`. In that file:

`overwrite` - `true` will replace the default properties with the external ones. `false` will append external properties to the end of the list.\
`endpoints` - list of url endpoints to fetch properties from.

#### `endpoints`

An endpoint is a url that can be queried for chemical properties. Each endpoint can be used to populate multiple properties. The currently supported modes are `GET smiles`, `POST smiles`, and `POST sdf`, with a response type of `json`.

In `GET smiles` mode, the endpoint `url` is expected to contain the string `:smiles`. When this endpoint is queried, the smiles for the chemical in question will replace `:smiles`. Example: `example.com/chem/:smiles` or `example.com/chem?q=:smiles`.

In `POST smiles` mode, the `payload` is expected to contain the string `:smiles`. When this endpoint is queried, the smiles for the chemical in question will replace `:smiles`. Example: `{"smiles":":smiles"}`.

`name` - unique name for the endpoint\
`url` - endpoint url\
`method` - either `GET` or `POST`\
`data` - either `sdf` or `smiles`\
`payload` - for `POST smiles` requests, the json payload to send containing `:smiles`\
`properties` - a mapping of property names to config

##### `properties`

The properties configuration allows one endpoint to populate multiple properties. Each property has a `path` to find the property in the response body. For example, if the endpoint response body is `{"properties":{"prop1":1}}`, the `path` would be `properties.prop1`.

`format` - format string to format the unit. example `"%.3f"`\
`path` - path to the property in the response body

## Development
To run Data Table with autoreload:

```sh
$ python3 -m pip install -r requirements.txt
$ python3 run.py -r [run_args]
```

## License

MIT
