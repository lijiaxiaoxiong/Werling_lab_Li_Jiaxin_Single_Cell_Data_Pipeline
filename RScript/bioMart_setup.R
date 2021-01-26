# Jan 18th
library(biomaRt)
listMarts()
Sys.setenv("http_proxy" = "http://my.proxy.org:9999")
options(RCurlOptions = list(proxy="uscache.kcc.com:80",proxyuserpwd="------:-------"))
