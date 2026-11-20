# QLever

NOTE: to run qlever as `root` user, some OS securities have to be disabled (`AppArmor`).
So better to use another user.

## QLever installation

### The client

```bash
uv tool install qlever
# or to update it
uv tool update qlever
```

### The container (tools)

```bash
podman pull docker.io/adfreiburg/qlever:latest  # same command to update
```

## Prepare Bgee data

Get ttl files from [https://www.bgee.org/ftp/16-0/rdf_easybgee.zip](https://www.bgee.org/ftp/16-0/rdf_easybgee.zip)

Bgee SIB SPARQL examples are available at [https://github.com/sib-swiss/sparql-examples](https://github.com/sib-swiss/sparql-examples)

```bash
unzip rdf_easybgee.zip
git clone https://github.com/sib-swiss/sparql-examples.git
cat rdf_data/*.ttl sparql-examples/examples/Bgee/[0-9]*.ttl >bgee.ttl
# or if too many files
find rdf_data/ -type f -name \*.ttl -exec cat {} \; > bgee.ttl; cat sparql-examples/examples/Bgee/[0-9]*.ttl >>bgee.ttl
gzip bgee.ttl
```

## Index ttl

To index `ttl` for QLever, you have to be in the folder where the `Qleverfile` file is located !!!

Also, qlever run a container, so it cannot see symlinks (e.g. the ttl file) linked outside the running directory !!!

```bash
cd $QLEVER_HOME
qlever stop
time qlever  index  --parallel-parsing false --overwrite-existing --stxxl-memory 50G
```

QLever options could be optimised on the command line AND in the `Qleverfile` file!

According to Jerven, the void data have to be loaded afterward: "I use the void generator at the end and then use a sparql update load to insert it".

## Run the endpoint

To run QLever, you have to be in the folder where the `Qleverfile` file is located !!!

```bash
cd $QLEVER_HOME
qlever start --access-token 'XYZ'
```

### The void generator

Download the jar file of the void generator from https://github.com/sib-swiss/void-generator

```bash
ln -s void-generator-X.YZ-uber.jar void-generator-latest.jar
java -jar void-generator-latest.jar  -f -r http://localhost:7030 --void-file=bgee_v16_0_void-generator.ttl --iri-of-void=https://www.bgee.org/sparql/.well-known/void
perl -i -pe "s|http://localhost:7030|https://www.bgee.org/sparql/|g" void.ttl # OK only valid for prod
```


### Daemonize/service

TODO

