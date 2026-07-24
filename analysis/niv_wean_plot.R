library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

# ── 模拟数据 ────────────────────────────────────────────────
day <- c(1, 3, 5, 7, 9, 11)

df <- data.frame(
  day    = day,
  PaO2   = c(55, 65, 76, 79, 81, 82),      # mmHg
  PaCO2  = c(82, 76, 72, 70, 69, 69),      # mmHg
  SpO2   = c(82, 86, 88, 89, 90, 90),      # %
  RR     = c(28, 24, 22, 20, 19, 18),      # 次/分
  HR     = c(112, 98, 92, 88, 84, 80),     # 次/分
  IPAP   = c(18, 16, 14, 12, 10, NA),      # cmH2O（撤机后NA）
  FiO2   = c(50, 45, 40, 35, 33, 28)       # %
)

# ── 脱机阶段标注 ────────────────────────────────────────────
phases <- data.frame(
  xmin  = c(1,   7,   9),
  xmax  = c(7,   9,  11),
  label = c("持续NIV治疗", "间断脱机\n（白天鼻导管）", "完全撤机\n鼻导管吸氧"),
  fill  = c("#FFF3CD", "#D1ECF1", "#D4EDDA")
)

phase_rect <- function(p) {
  p +
    annotate("rect", xmin=1, xmax=7,  ymin=-Inf, ymax=Inf, alpha=0.15, fill="#FFC107") +
    annotate("rect", xmin=7, xmax=9,  ymin=-Inf, ymax=Inf, alpha=0.15, fill="#17A2B8") +
    annotate("rect", xmin=9, xmax=11, ymin=-Inf, ymax=Inf, alpha=0.15, fill="#28A745") +
    geom_vline(xintercept=c(7,9), linetype="dashed", color="gray50", linewidth=0.4)
}

base_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size=10, face="bold"),
    axis.title.x = element_blank()
  )

x_scale <- scale_x_continuous(breaks=day, labels=paste0("D", day))

# ── 图1：血气 PaO2 & PaCO2 ──────────────────────────────────
p1 <- ggplot(df, aes(x=day)) +
  geom_line(aes(y=PaO2,  color="PaO₂"),  linewidth=1.2) +
  geom_point(aes(y=PaO2, color="PaO₂"),  size=3) +
  geom_line(aes(y=PaCO2, color="PaCO₂"), linewidth=1.2) +
  geom_point(aes(y=PaCO2,color="PaCO₂"), size=3) +
  geom_hline(yintercept=60, linetype="dotted", color="#E74C3C", linewidth=0.5) +
  annotate("text", x=11.2, y=60, label="PaO₂目标≥60", size=3, color="#E74C3C", hjust=1) +
  scale_color_manual(values=c("PaO₂"="#2196F3","PaCO₂"="#FF5722")) +
  labs(title="① 动脉血气（mmHg）", y="mmHg", color=NULL) +
  x_scale + base_theme +
  theme(legend.position="top")
p1 <- phase_rect(p1)

# ── 图2：SpO2 ────────────────────────────────────────────────
p2 <- ggplot(df, aes(x=day, y=SpO2)) +
  geom_line(color="#9C27B0", linewidth=1.2) +
  geom_point(color="#9C27B0", size=3) +
  geom_hline(yintercept=88, linetype="dotted", color="#E74C3C", linewidth=0.5) +
  annotate("text", x=11.2, y=88, label="目标≥88%", size=3, color="#E74C3C", hjust=1) +
  scale_y_continuous(limits=c(80,93), breaks=seq(80,93,2)) +
  labs(title="② 指氧饱和度 SpO₂（%）", y="%") +
  x_scale + base_theme
p2 <- phase_rect(p2)

# ── 图3：呼吸频率 & 心率 ────────────────────────────────────
df_vital <- df %>% select(day, RR, HR) %>%
  pivot_longer(c(RR,HR), names_to="param", values_to="val")

p3 <- ggplot(df_vital, aes(x=day, y=val, color=param)) +
  geom_line(linewidth=1.2) +
  geom_point(size=3) +
  geom_hline(yintercept=20, linetype="dotted", color="#E74C3C", linewidth=0.5) +
  annotate("text", x=11.2, y=20, label="RR目标≤20", size=3, color="#E74C3C", hjust=1) +
  scale_color_manual(values=c("RR"="#FF9800","HR"="#607D8B"),
                     labels=c("RR 呼吸频率（次/分）","HR 心率（次/分）")) +
  labs(title="③ 生命体征", y="次/分", color=NULL) +
  x_scale + base_theme +
  theme(legend.position="top")
p3 <- phase_rect(p3)

# ── 图4：呼吸机参数 IPAP & FiO2 ─────────────────────────────
p4 <- ggplot(df, aes(x=day)) +
  geom_line(aes(y=IPAP, color="IPAP"), linewidth=1.2, na.rm=TRUE) +
  geom_point(aes(y=IPAP,color="IPAP"), size=3, na.rm=TRUE) +
  geom_point(aes(x=11, y=NA), na.rm=TRUE) +  # 撤机点
  annotate("point", x=11, y=6, shape=4, size=5, color="#4CAF50", stroke=2) +
  annotate("text",  x=11, y=7, label="撤机", size=3.5, color="#4CAF50") +
  geom_line(aes(y=FiO2, color="FiO₂"), linewidth=1.2) +
  geom_point(aes(y=FiO2,color="FiO₂"), size=3) +
  scale_color_manual(values=c("IPAP"="#3F51B5","FiO₂"="#009688")) +
  labs(title="④ 呼吸机参数", y="cmH₂O / %", color=NULL,
       caption="IPAP单位：cmH₂O；FiO₂单位：%；Day11 IPAP=撤机") +
  x_scale + base_theme +
  theme(legend.position="top",
        plot.caption=element_text(size=8, color="gray50"))
p4 <- phase_rect(p4)

# ── 拼图 ─────────────────────────────────────────────────────
combined <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title    = "AECOPD患者无创呼吸机脱机过程监测趋势",
    subtitle = "黄色：持续NIV阶段 | 蓝色：间断脱机阶段 | 绿色：完全撤机阶段",
    theme = theme(
      plot.title    = element_text(size=13, face="bold"),
      plot.subtitle = element_text(size=10, color="gray40")
    )
  )

ggsave("analysis/figures/NIV_wean_monitoring.png",
       combined, width=13, height=9, dpi=180)

message("图已保存至 analysis/figures/NIV_wean_monitoring.png")

# ── 打印模拟数据表 ───────────────────────────────────────────
cat("\n=== 模拟监测数据 ===\n")
df_print <- df
df_print$IPAP[is.na(df_print$IPAP)] <- "撤机"
df_print$day <- paste0("Day", df_print$day)
print(df_print, row.names=FALSE)
