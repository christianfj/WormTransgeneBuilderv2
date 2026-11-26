# Running the app as a service

The app is distributed as a Docker image. Any server that can run a current Docker can host it. These are the
steps used for the public site (2026-11-26); host names and addresses are left out on purpose.

## Requirements

* A server with a **current Linux** (Ubuntu 24.04 or later) and a **current Docker** (Docker 20.10 or later; the version in
  Ubuntu's own package archive is fine). Older systems cannot run the image: Ubuntu 16.04 with Docker 18.09 stops it with the
  error `R_HOME ('/usr/local/lib/R') not found`, because that Docker's security filter does not know a system call that modern
  R uses. Do not "fix" this by switching the filter off (`--security-opt seccomp=unconfined`); update Docker instead.
* 2 CPUs and 4 GB of memory is a minimum; 8 GB or more is comfortable. A 20 kb design is the heavy case (the piRNA comparison uses all
  CPUs for several seconds, and the GLO scores take about 130 MB of memory per R process).
* About 5 GB of disk for Docker and the image (1.7 GB).
* Ports 80 and 443 open to the internet (the web server); everything else closed (see Firewall).

## Start the app

```sh
docker run -d --name wormtransgenebuilder --restart unless-stopped --memory 6g \
  -p 127.0.0.1:3839:3838 ghcr.io/christianfj/wormtransgenebuilderv2:v2.0.0
```

* `--restart unless-stopped` starts it again after a reboot or a crash.
* `--memory` stops one heavy request from using the whole machine. 6 GB was chosen for a 16 GB machine that also runs other
  services; adjust it to your machine.
* `-p 127.0.0.1:3839:3838` publishes the port **on the server only**. The app is reached through the web server below and
  never directly from the internet. The container's port is always 3838; 3839 is any free port on the host.
* The image is public: no login is needed to pull it. The container runs as an unprivileged user and has a health check
  (`docker ps` shows `healthy`). Its log starts with a line giving the app, R and ViennaRNA versions
  (`docker logs wormtransgenebuilder`), and the **About** tab shows the same.

Check it from the server: `curl -sI http://127.0.0.1:3839/` should answer `200 OK`.

## Put it on the web with HTTPS (nginx)

The public address is a **sub-path** (`/transgenebuilder/`) of an existing site that already has HTTPS and a certificate. Add this
inside that site's `server { ... }` block:

```
location = /transgenebuilder { return 301 /transgenebuilder/; }
location /transgenebuilder/ {
  proxy_pass http://127.0.0.1:3839/;        # the trailing slash removes the /transgenebuilder/ prefix
  proxy_http_version 1.1;
  proxy_set_header Upgrade $http_upgrade;   # Shiny needs WebSockets
  proxy_set_header Connection $connection_upgrade;
  proxy_set_header Host $host;
  proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
  proxy_set_header X-Forwarded-Proto $scheme;
  proxy_read_timeout 20d;                   # an open session is a long-lived WebSocket
}
```

`$connection_upgrade` must be defined once, at the top level of the nginx configuration (outside any `server` block):

```
map $http_upgrade $connection_upgrade {
  default upgrade;
  ''      close;
}
```

To serve the app at the root of its own domain instead, use `location / { proxy_pass http://127.0.0.1:3839; ... }` with the same header lines
and let a tool such as Caddy or certbot handle the certificate.

Safe way to change it: back up the file, edit, then `sudo nginx -t` (it must say the test is successful), then
`sudo systemctl reload nginx` (a reload does not drop running connections; do not restart). To go back, restore the backup and reload.
Test with `curl -sI https://your.domain/transgenebuilder/` and by running a design in a browser: a page that stays grey or loses its
connection means the WebSocket lines or the trailing slash are wrong.

## Firewall

Only ports 80 and 443 should be reachable from the whole internet. Allow SSH (port 22) only from the addresses you administer
from. Close any port you used for testing (the app itself needs none, because it is bound to the server's own address).

## Updating

```sh
docker pull ghcr.io/christianfj/wormtransgenebuilderv2:<new version>
docker rename wormtransgenebuilder wormtransgenebuilder-previous
docker stop wormtransgenebuilder-previous
docker run -d --name wormtransgenebuilder --restart unless-stopped --memory 6g \
  -p 127.0.0.1:3839:3838 ghcr.io/christianfj/wormtransgenebuilderv2:<new version>
# check it, then remove the old container; if something is wrong, stop the new one and "docker start wormtransgenebuilder-previous"
docker rm wormtransgenebuilder-previous
```

Keep the previous image until the new one has been checked, so you can go back. Enable automatic security updates of the server's
operating system (`unattended-upgrades` on Ubuntu) and restart the server when it asks for it (the container comes back by itself).

## Things to know about cloud servers

* On some Amazon EC2 types (for example the `i3`, `i4i` and `d` families) a large local SSD is attached. It is **instance store**:
  it is wiped when the instance is *stopped* or terminated, or if the hardware fails (a reboot keeps it). If Docker's data lives there,
  the image has to be pulled again after a stop, which is harmless because the image is public, but keep nothing irreplaceable on it.
  A plain `t3` or similar instance with a normal disk has no such trap.
* Give the server a fixed address (an Elastic IP on AWS) if the domain name points to it.
* Set a billing alert so that a surprise shows up before the bill does.

## Without a registry

`docker save` / `docker load` (see [REPRODUCIBILITY.md](REPRODUCIBILITY.md)) moves an image between machines
without any registry. Every release carries a saved copy of its image, and the same copy is archived on Zenodo.
