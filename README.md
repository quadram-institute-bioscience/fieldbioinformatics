### Coronahit
This branch `coronahit` is a small modification of the artic pipeline [1] in which we use `samtools ampliconclip` for trimming primers instead of `align_trim` coming along with the pipeline. As the subcommand `ampliconclip` is only available in samtools >=1.11 and some dependencies for artic pipeline require older samtools version, so we currently use docker with a compiled samtools 1.11 and a conda env for `artic` as a workaround. Therefore, to use the modified artic for coronahit protocol, the most convenient way is to run it via docker

Assuming you have data in a directory in a Linux machine `/home/ubuntu/coronahit`, the command for running artic: 

```
docker run --rm -v /home/ubuntu/coronahit:/data quadram/artic-coronahit artic --help
```
`/data` is a predefined working directory inside the docker image `quadram/artic-coronahit`

**Notes**

If you use Docker Desktop on MacOS or Windows, make sure you allocate suitable computing resources for docker. See [this manual](https://docs.docker.com/docker-for-windows/#resources) for example.

[1] https://github.com/artic-network/fieldbioinformatics