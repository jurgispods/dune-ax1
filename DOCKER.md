# Run simulations in Docker image

Compiling an Autotools-managed C++ application with multiple dependencies and _very_ specific version requirements is hard. So I created a Docker image specification for ease-of-use and portability.

## Prerequites

You need 

- a somewhat capable computer
- and a Docker (and optionally, `docker-compose`) installation 
which has been granted around 8-16GB of RAM and as many vCPUs as you can spare :)

You have the following options:
- Pull the official stable `dune-ax1` Docker image und run without modification
- Pull a 'development' image which provides the complete build environment and therefore allows code modifications and recompilation from within the image
- Build the images yourself

### Official Docker image

```
docker pull pederpansen/dune-ax1
```

### Dev Docker image

```
docker pull pederpansen/dune-ax1:dev
```

### Build images

Build dev image:

```
docker build -t pederpansen/dune-ax1:dev -f Dockerfile_build .
```

Build stable image:

```
docker build -t pederpansen/dune-ax1:latest .
```

## Run simulation

```
docker-compose run --rm dune-ax1
```

This will run the setup from my [Biophysical Journal 2013 paper](https://doi.org/10.1016/j.bpj.2013.05.041)

## Copy simulations results from container to host


```
docker cp <container name>:/dune/output .
```

This is necessary mostly when running under Mac or Windows for performance reasons, where the Docker container does
not write to the host filesystem, but to a VM filesystem.

If you are running under Linux, you might as well uncomment the relevant line in `docker-compose.yml` and directly bind-mount a host folder to the container.
