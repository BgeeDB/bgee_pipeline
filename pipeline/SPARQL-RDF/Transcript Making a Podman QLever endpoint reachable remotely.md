# Transcript: Making a Podman QLever endpoint reachable remotely

Nov 27, 2026 · @Sébastien Moretti

**Summary:** QLever runs in a root Podman container on port 7030. It could not be reached from another machine, while Virtuoso on port 9940 on the same server could. The cause was UFW's default DROP policy for forwarded (routed) traffic: `ufw allow 7030` only opens the port for processes on the host, not for traffic forwarded into a container. The fix was two `ufw route allow` rules.

---

## 1. User

How can I test a local SPARQL endpoint on the command line?

## 1. Claude

The quickest way is `curl`, since SPARQL over HTTP is just a GET or POST request. Swap in your own endpoint URL; common defaults are `http://localhost:3030/<dataset>/sparql` for Fuseki, `http://localhost:8890/sparql` for Virtuoso, `http://localhost:7200/repositories/<repo>` for GraphDB, and `http://localhost:9999/blazegraph/sparql` for Blazegraph.

**1. Basic query (POST, form-encoded)**

```bash
curl -s http://localhost:3030/ds/sparql \
  -H 'Accept: application/sparql-results+json' \
  --data-urlencode 'query=SELECT * WHERE { ?s ?p ?o } LIMIT 5'
```

**2. GET request**

```bash
curl -sG http://localhost:3030/ds/sparql \
  -H 'Accept: application/sparql-results+json' \
  --data-urlencode 'query=SELECT * WHERE { ?s ?p ?o } LIMIT 5'
```

The `-G` flag makes curl put the URL-encoded data into the query string.

**3. Send a query from a file**

```bash
curl -s http://localhost:3030/ds/sparql \
  -H 'Content-Type: application/sparql-query' \
  -H 'Accept: text/csv' \
  --data-binary @query.rq
```

**4. Simple health check**

```bash
curl -s -o /dev/null -w '%{http_code}\n' \
  --data-urlencode 'query=ASK { ?s ?p ?o }' \
  http://localhost:3030/ds/sparql
```

A `200` means the endpoint is up. Drop `-o /dev/null` and you'll see `true` if the store has any data.

**5. Updates** (these usually go to a separate `/update` path)

```bash
curl -s http://localhost:3030/ds/update \
  --data-urlencode 'update=INSERT DATA { <urn:a> <urn:b> "test" }'
```

**Choosing the output format** with the `Accept` header:

- `application/sparql-results+json` for JSON. Pipe it to `jq` to make it readable: `| jq '.results.bindings'`
- `text/csv` or `text/tab-separated-values` for tables
- `text/turtle` for CONSTRUCT and DESCRIBE queries

**Dedicated tools**, if you'd rather not use curl:

- **Apache Jena** includes `rsparql --service http://localhost:3030/ds/sparql --query query.rq --results text`, which prints a nicely aligned table.
- **Python:** run `pip install sparqlwrapper` for scripted tests.
- **Comunica:** run `npx @comunica/query-sparql http://localhost:3030/ds/sparql "SELECT * WHERE {?s ?p ?o} LIMIT 5"`

**Troubleshooting:**

- **Connection refused:** the server isn't running, or it's on a different port.
- **404:** the dataset or repository name in the path is wrong.
- **400:** usually a syntax error in the query, or you used `-d` instead of `--data-urlencode`, so characters like `{` or `#` weren't encoded.
- **401 or 403:** add `-u user:pass`.

---

## 2. User

On the same machine, I have virtuoso on port 9940 and QLever on port 7030. On another machine, I can query virtuoso but not QLever. The firewall on the first machine is opened for both ports.

## 2. Claude

Virtuoso works from the other machine, so the network and firewall are fine in general. The difference is in how QLever is set up. Work through these checks in order.

**1. On the remote machine, look at the error you get**

```bash
curl -v 'http://SERVER:7030/?query=ASK%7B%3Fs%3Fp%3Fo%7D'
```

- **`Connection refused`:** nothing is listening on that address. QLever is probably bound to localhost only.
- **Timeout, with no reply at all:** something is dropping the packets. This is usually firewall or Docker iptables rules.
- **An HTTP error:** QLever is reachable, and the problem is the request itself.

**2. On the server, check which address QLever listens on**

```bash
ss -tlnp | grep -E '7030|9940'
```

You want `0.0.0.0:7030`, `*:7030` or `[::]:7030`. If it says `127.0.0.1:7030`, it only accepts local connections. If QLever runs in Docker, also run `docker ps --format '{{.Names}}\t{{.Ports}}'`. If you see `127.0.0.1:7030->7030/tcp`, re-run the container with `-p 7030:7030`.

**3. Docker and the firewall**

Docker manages its own iptables chains, so "port opened in the firewall" doesn't always apply to containers. With firewalld, Docker traffic can go through the `docker` zone instead of your public zone; look for DROP rules in `iptables -L DOCKER-USER -n` (or `nft list ruleset`). With rootless Podman, port forwarding works differently. A quick test is to curl the server's LAN IP from the server itself: if that fails but `localhost:7030` works, the cause is binding or Docker, not the firewall.

**4. Request differences between the two engines**

- **Path:** QLever serves SPARQL at the root (`http://host:7030/`), while Virtuoso uses `/sparql`.
- **Updates:** QLever only accepts updates with an access token.
- **Proxy:** if the remote machine has `http_proxy` set, the proxy may allow 9940 but not 7030. Test with `curl --noproxy '*' ...`.

---

## 3. User

```
$ curl -v 'http://130.223.29.25:7030/?query=ASK%7B%3Fs%3Fp%3Fo%7D'
*   Trying 130.223.29.25:7030...
* connect to 130.223.29.25 port 7030 from 130.223.28.162 port 45860 failed: Connection timed out
* Failed to connect to 130.223.29.25 port 7030 after 134599 ms: Could not connect to server
curl: (28) Failed to connect to 130.223.29.25 port 7030 after 134599 ms: Could not connect to server

$ ss -tlnp | grep -E '7030|9940'
LISTEN 0      500          0.0.0.0:9940       0.0.0.0:*
LISTEN 0      4096         0.0.0.0:7030       0.0.0.0:*

# podman ps --format '{{.Names}}\t{{.Ports}}'
qlever.server.bgee  0.0.0.0:7030->7030/tcp
qlever.ui.bgee      0.0.0.0:7031->7000/tcp
```

## 3. Claude

The binding is fine: QLever listens on `0.0.0.0:7030`, and Podman publishes it on all interfaces. The **timeout** means the connection packets are being silently dropped.

**Most likely cause: the firewall was reloaded after the containers started.** Virtuoso runs natively, so incoming traffic only passes the INPUT rules. Traffic to a Podman container is DNAT'ed and *forwarded*, using rules that Podman adds when the container starts. A firewall reload wipes them. Fix, as root:

```bash
podman network reload --all
```

**If that doesn't fix it, check whether the packets reach the server at all:**

```bash
tcpdump -ni any 'tcp port 7030'
```

- **No packets appear:** they're blocked before reaching the server, e.g. by a network ACL between the `130.223.28.x` and `130.223.29.x` subnets.
- **SYN packets arrive but nothing goes to the container interface:** the host is dropping forwarded traffic. Check `firewall-cmd --get-active-zones`, `firewall-cmd --list-all` and `nft list ruleset | grep -n 7030`.

---

## 4. User

```
# podman network reload --all
a316aefcb6f9a4671a5c4cd476536c9d4b726689aa3bbf1a471f763ea447bd09
fc1595c0e90f53d5ae58f26c96bf3979d1afbd60244c6653e224bcf8bbb94d0d

$ curl -v 'http://130.223.29.25:7030/?query=ASK%7B%3Fs%3Fp%3Fo%7D'
* connect to 130.223.29.25 port 7030 from 130.223.28.162 port 55732 failed: Connection timed out
curl: (28) Failed to connect to 130.223.29.25 port 7030 after 132762 ms: Could not connect to server

# tcpdump -ni any 'tcp port 7030'
17:06:10.061509 ens801f0np0 In  IP 130.223.28.162.55732 > 130.223.29.25.7030: Flags [S], ...
17:06:10.061509 bond0 In  IP 130.223.28.162.55732 > 130.223.29.25.7030: Flags [S], ...
17:06:22.490683 ens801f1np1 In  IP 10.67.20.27.57710 > 130.223.29.25.7030: Flags [S], ...
17:06:22.490683 bond0 In  IP 10.67.20.27.57710 > 130.223.29.25.7030: Flags [S], ...
(SYN retransmissions only, no replies)

# nft list ruleset | grep -n 7030
168:		tcp dport 7030 counter packets 0 bytes 0 accept
169:		udp dport 7030 counter packets 0 bytes 0 accept
410:		tcp dport 7030 counter packets 0 bytes 0 accept
411:		udp dport 7030 counter packets 0 bytes 0 accept
456:		ip protocol tcp  tcp dport 7030 counter packets 22 bytes 1348 jump CNI-DN-0f9006a2534928135f678
486:		ip saddr 10.88.0.0/16 tcp dport 7030 counter packets 0 bytes 0 jump CNI-HOSTPORT-SETMARK
487:		ip saddr 127.0.0.1 tcp dport 7030 counter packets 0 bytes 0 jump CNI-HOSTPORT-SETMARK
488:		tcp dport 7030 counter packets 22 bytes 1348 dnat to 10.88.0.14:7030
```

## 4. Claude

The server receives the SYN packets, and the NAT rule rewrites them to the container (`dnat to 10.88.0.14:7030`, 22 packets). After that they disappear, so the server drops them in the **forwarding** step. The port you opened (lines 168–169) has 0 packets because it's in the INPUT chain; after the address rewrite, container traffic goes through the **FORWARD** chain.

**1. Check that IP forwarding is enabled:** `sysctl net.ipv4.ip_forward`. If `0`, enable it and persist it in `/etc/sysctl.d/`.

**2. If forwarding is already `1`, look at the FORWARD chain:** `iptables -L FORWARD -n -v --line-numbers`. The likely culprit is a DROP policy or rule before the jump to `CNI-FORWARD`. Allow forwarding to and from the Podman bridge (usually `cni-podman0`, subnet `10.88.0.0/16`), for ports 7030 and 7000 (the UI's internal port).

**Alternative:** start the QLever containers with `--network host` so QLever binds directly to the host's port, like Virtuoso.

The capture also shows a second client, `10.67.20.27`, failing the same way; the fix should solve it too.

---

## 5. User

```
# sysctl net.ipv4.ip_forward
net.ipv4.ip_forward = 1

# iptables -L FORWARD -n -v --line-numbers
Chain FORWARD (policy DROP 289 packets, 17816 bytes)
num   pkts bytes target     prot opt in     out     source               destination
1      289 17816 CNI-FORWARD  ...  /* CNI firewall plugin rules */
2      289 17816 ufw-before-logging-forward  ...
3      289 17816 ufw-before-forward  ...
4      289 17816 ufw-after-forward  ...
5      289 17816 ufw-after-logging-forward  ...
6      289 17816 ufw-reject-forward  ...
7      289 17816 ufw-track-forward  ...
```

Also, QLever runs as root, with a service/daemon.

## 5. Claude

The firewall on this server is **UFW**, and its default policy for forwarded ("routed") traffic is DROP. All 289 packets went through every chain without being accepted, then hit the `policy DROP`.

`ufw allow 7030` only opens the port for processes running on the host, like Virtuoso. Traffic forwarded into a Podman container needs a **route** rule:

```bash
ufw route allow proto tcp from any to 10.88.0.0/16 port 7030
ufw route allow proto tcp from any to 10.88.0.0/16 port 7000
```

The second rule is for the QLever UI (host port 7031 → container port 7000). Route rules are checked after the address is rewritten, so they use the container ports and the Podman subnet; using the subnet keeps them working when the container gets a new IP. Retest with:

```bash
curl -s 'http://130.223.29.25:7030/?query=ASK%7B%3Fs%3Fp%3Fo%7D'
```

UFW saves route rules itself, so they survive reboots and service restarts (check with `ufw status numbered`, shown as `ALLOW FWD`). Setting `DEFAULT_FORWARD_POLICY="ACCEPT"` in `/etc/default/ufw` would also work but opens every container port. Running QLever as a root service doesn't change anything here.

---

## 6. User

Thanks It works! I need a transcript of this chat.
