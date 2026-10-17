wdiag <- function(ccw, version) {
    caps <- ccw %>% group_by(strategy) %>% summarise(cap = unname(quantile(ipcw_w, 
        0.99, na.rm = TRUE)), .groups = "drop")
    d <- ccw %>% left_join(caps, by = "strategy")
    d %>% group_by(strategy) %>% summarise(version = version, 
        cap_p99_full_arm = first(cap), trimmed_pct_full_arm = 100 * 
            mean(ipcw_w > cap), n_equipoise = sum(equipoise %in% 
            TRUE), trimmed_pct_equipoise = 100 * mean(ipcw_w[equipoise %in% 
            TRUE] > cap[equipoise %in% TRUE]), ess_equipoise = {
            w <- ipcw_w_trim[equipoise %in% TRUE]
            sum(w)^2/sum(w^2)
        }, mean_w_equipoise = mean(ipcw_w_trim[equipoise %in% 
            TRUE]), max_w_equipoise = max(ipcw_w_trim[equipoise %in% 
            TRUE]), .groups = "drop")
}
