# Running the app as a service

The app is distributed as a Docker image. Any server that can run Docker can host it.

```sh
docker run -d --name wormbuilder --restart unless-stopped \
  --memory 4g --cpus 2 -p 127.0.0.1:3838:3838 \
  ghcr.io/christianfj/wormtransgenebuilderv2:v2.0.0
```

* `--restart unless-stopped` starts it again after a reboot or a crash.
* `--memory` and `--cpus` stop one heavy request from using the whole machine (a 20 kb design is the heavy case; GLO scoring keeps about 130 MB in memory per R process).
* Publishing the port on `127.0.0.1` only means the app is reached through a reverse proxy.
* The container runs as an unprivileged user and has a health check (`docker ps` shows `healthy`).

## HTTPS

Put a reverse proxy in front for HTTPS and certificates (for example Caddy or nginx with Let's Encrypt).
Shiny uses WebSockets, so the proxy must allow them. Example with Caddy:

```
your.domain.example {
    reverse_proxy 127.0.0.1:3838
}
```

## Updating

```sh
docker pull ghcr.io/christianfj/wormtransgenebuilderv2:<new version>
docker rm -f wormbuilder
docker run ... (as above, with the new tag)
```

Keep the previous tag so you can go back. Enable automatic security updates of the server's operating
system (`unattended-upgrades` on Ubuntu).

## Without a registry

`docker save` / `docker load` (see [REPRODUCIBILITY.md](REPRODUCIBILITY.md)) moves an image between machines
without any registry.
