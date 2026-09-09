# -*- coding: utf-8 -*-
"""
ARDS Bedside Physiological Calculator & Dead-Space Estimator v11
臨床床邊生理指標、死腔預估、可復張性(R/I Ratio)、自主呼吸驅力與生物表型預測計算器 (Streamlit Web App & CLI 雙模工具)
Based on:
1. Nuckton 2002 (NEJM) & Sinha 2019 (AJRCCM) - Ventilatory Ratio (VR) & Dead Space
2. Wongtirawit 2026 (ICM) & Chen 2020 (AJRCCM) - R/I Ratio, Ventilation Intensity, AOP Correction
3. Sinha 2020 (Lancet Respiratory Medicine) - Biological Phenotype Parsimonious Models
"""

import sys
import math

# ---------------------------------------------------------
# Physiological Core Calculations
# ---------------------------------------------------------

def calculate_pbw(height_cm, sex):
    """Calculate Ideal Body Weight / Predicted Body Weight (PBW) in kg."""
    if sex.lower() in ['male', 'm', '男'] :
        return 50.0 + 0.91 * (height_cm - 152.4)
    else:
        return 45.5 + 0.91 * (height_cm - 152.4)

def calculate_vr(ve, paco2, pbw):
    """Calculate Ventilatory Ratio (VR)."""
    if pbw <= 0:
        return 0
    return (ve * paco2) / (4.0 * pbw)

def estimate_vd_vt(vr):
    """Estimate Physiological Dead-Space Fraction (Vd/Vt) from VR."""
    if vr <= 0.7:
        return 0.30  # Floor value for normal physiological dead space
    return 1.0 - (0.70 / vr)

def get_risk_status(vr, vd_vt):
    """Determine risk category based on Nuckton 2002 quintile thresholds."""
    if vr < 1.50:
        return {
            "level": "GREEN",
            "title": "低風險區 (Low Risk)",
            "color": "#2ebd59",
            "desc": "通氣效率良好，死腔比例低於 ARDS 存活者平均值 (Vd/Vt < 0.54)。",
            "action": "維持當前保護性通氣策略。持續監測呼吸力學。"
        }
    elif vr < 1.90:
        return {
            "level": "YELLOW",
            "title": "中度警戒區 (Intermediate Risk)",
            "color": "#f1c40f",
            "desc": "死腔比例介於存活組與死亡組平均值之間 (Vd/Vt 0.54 ~ 0.63)。死亡風險隨死腔比例增加而上升。",
            "action": "1. 密切監測 PEEP 調整對 VR 的動態影響，尋找最佳 PEEP。\\n2. 注意有無小氣道塌陷或早期肺部順應性變差。"
        }
    else:
        # Severe or extreme
        if vr >= 2.25:
            return {
                "level": "RED_EXTREME",
                "title": "極高度危殆區 (Extreme Risk - Quintile 5)",
                "color": "#9b0000",
                "desc": "死腔比例極高 (Vd/Vt >= 0.69，落入 Nuckton 研究最嚴重的第五個五分位區間)。",
                "action": "1. 強烈建議排除肺微血管微血栓 (Microthrombi) 或血管病變。\\n2. 立即安排心臟超音波評估右心功能，防範急性肺心症 (ACP)。\\n3. 考慮實施更嚴格的肺與右心保護 (Pplat < 26-28 cmH2O, pH > 7.25)。"
            }
        else:
            return {
                "level": "RED",
                "title": "高度危險區 (High Risk)",
                "color": "#e74c3c",
                "desc": "死腔比例已達或超過 ARDS 死亡組 the average value (Vd/Vt >= 0.63)。通氣效率嚴重受損。",
                "action": "1. 評估右心負荷與肺血管阻力 (PVR)。\\n2. 優化水分平衡管理 (FACCT 限制性補水/利尿)。\\n3. 評估是否需啟動早期俯臥通氣 (Prone) 以優化 V/Q 匹配。"
            }

# 🔄 R/I Ratio Calculations
def calculate_ri_ratio(peep_high, peep_low, vt_high, v_exp_drop, pplat_low, vt_low):
    """
    Calculate Recruitment-to-Inflation (R/I) ratio.
    """
    if (pplat_low - peep_low) <= 0:
        return 0.0, 0.0, 0.0, 0.0
    c_low = vt_low / (pplat_low - peep_low)
    delta_peep = peep_high - peep_low
    delta_eelv = v_exp_drop - vt_high
    v_recruited = delta_eelv - (c_low * delta_peep)
    if (c_low * delta_peep) <= 0:
        return 0.0, c_low, delta_eelv, 0.0
    ri_ratio = v_recruited / (c_low * delta_peep)
    return ri_ratio, c_low, delta_eelv, v_recruited

# 🏆 Primary Phenotype Models (Figure 2, Sinha 2020)
def calculate_primary_il8_4var_prob(il8, protein_c, bicarbonate, vasopressor_used):
    vp = 1.0 if vasopressor_used else 0.0
    log_il8 = math.log(il8 + 1.0)
    log_prot_c = math.log(protein_c + 1.0)
    y = 4.8127 + 1.5424 * log_il8 - 0.2502 * bicarbonate - 1.8778 * log_prot_c + 2.1026 * vp
    return 1.0 / (1.0 + math.exp(-y))

def calculate_primary_il8_3var_prob(il8, protein_c, bicarbonate):
    log_il8 = math.log(il8 + 1.0)
    log_prot_c = math.log(protein_c + 1.0)
    y = 6.1241 + 1.4226 * log_il8 - 0.2596 * bicarbonate - 1.8330 * log_prot_c
    return 1.0 / (1.0 + math.exp(-y))

# 🧬 Ancillary Phenotype Models (Supplementary Table S4)
def calculate_ancillary_il6_4var_prob(il6, protein_c, bicarbonate, vasopressor_used):
    vp = 1.0 if vasopressor_used else 0.0
    log_il6 = math.log(il6 + 1.0)
    log_prot_c = math.log(protein_c + 1.0)
    y = 4.0323 - 0.2581 * bicarbonate + 1.7412 * vp - 1.3805 * log_prot_c + 0.9191 * log_il6
    return 1.0 / (1.0 + math.exp(-y))

def calculate_ancillary_il8_sTNFR1_4var_prob(il8, stnfr1, bicarbonate, vasopressor_used):
    vp = 1.0 if vasopressor_used else 0.0
    log_il8 = math.log(il8 + 1.0)
    log_stnfr1 = math.log(stnfr1 + 1.0)
    y = -13.1351 + 1.3947 * log_il8 - 0.2145 * bicarbonate + 1.1818 * log_stnfr1 + 2.1398 * vp
    return 1.0 / (1.0 + math.exp(-y))

# ---------------------------------------------------------
# Streamlit App Mode
# ---------------------------------------------------------
def run_streamlit():
    import streamlit as st
    
    st.set_page_config(
        page_title="ARDS Bedside Physiological & Biological Phenotype Master Calculator",
        page_icon="🫁",
        layout="centered"
    )
    
    # Custom CSS for clinical layout
    st.markdown("""
    <style>
    .main-title {
        font-size: 26px;
        font-weight: bold;
        color: #2c3e50;
        text-align: center;
        margin-bottom: 5px;
    }
    .subtitle {
        font-size: 14px;
        color: #7f8c8d;
        text-align: center;
        margin-bottom: 25px;
    }
    .metric-container {
        background-color: #f8f9fa;
        padding: 15px;
        border-radius: 10px;
        border: 1px solid #e9ecef;
        text-align: center;
        margin-bottom: 10px;
    }
    .metric-value {
        font-size: 30px;
        font-weight: bold;
        color: #2c3e50;
    }
    .metric-label {
        font-size: 13px;
        color: #7f8c8d;
        font-weight: 500;
    }
    .warning-box {
        padding: 20px;
        border-radius: 10px;
        color: white;
        margin-top: 15px;
        margin-bottom: 15px;
    }
    .tab-content {
        padding-top: 15px;
    }
    .formula-box {
        background-color: #f1f2f6;
        padding: 12px;
        border-radius: 8px;
        border-left: 5px solid #3498db;
        font-family: monospace;
        margin-bottom: 15px;
        font-size: 13px;
        color: #2c3e50;
    }
    </style>
    """, unsafe_allow_html=True)
    
    st.markdown('<div class="main-title">🫁 ARDS 精準生理監測與生物亞型大師級計算器</div>', unsafe_allow_html=True)
    st.markdown('<div class="subtitle">整合 R/I Ratio、AOP 校正驅動壓、通氣強度、自主呼吸 efforts 與發炎表型評估</div>', unsafe_allow_html=True)
    
    # 5 Tabs Setup - Sinha Phenotype moved to Tab 1 (immediately after Tab 0)
    tab0, tab1, tab2, tab3, tab4 = st.tabs([
        "ℹ️ ARDS 臨床簡介與最新指引",
        "🧬 1. ARDS 生物亞型預測 (Sinha 2020)",
        "🫁 2. 呼吸力學與通氣效率", 
        "🔄 3. 肺可復張性評估 (R/I Ratio)", 
        "🧠 4. 自主呼吸 efforts 與驅力監測"
    ])
    

    # ---------------------------------------------------------
    # TAB 0: ARDS Clinical Introduction & Guidelines
    # ---------------------------------------------------------
    with tab0:
        st.markdown('<div class="tab-content"></div>', unsafe_allow_html=True)
        st.subheader("ℹ️ ARDS 臨床簡介與最新指引 (2024 Global ARDS Definition)")
        
        # Section A: 2024 Global Definition
        st.markdown("""
        ARDS 最新臨床定義於 **2024 年發布之「全球新定義 (A New Global Definition of ARDS)」**。
        此定義由 32 位重症醫學專家共同提出，發表於《American Journal of Respiratory and Critical Care Medicine》，是在 2012 年 Berlin 定義基礎上的重大擴展。
        其主要旨在解決臨床高流量鼻導管氧療（HFNO）普及、脈搏血氧飽和度（SpO2）取代動脈血氣的趨勢、超音波影像的應用，以及資源有限地區的床邊實用診斷需求。
        """)
        
        st.markdown(r"""
        ### 📌 2024 全球定義的四項核心建議：
        1. **納入高流量鼻導管氧療 (HFNO)**：
           - 只要 HFNO 流速 $\ge$ 30 L/min，即可在未插管狀態下直接診斷 ARDS，不再強求一定要 PEEP 或正壓通氣。
        2. **低血氧判定標準雙軌並行**：
           - 在血氧飽和度 $SpO_2 \le 97\%$ 的前提下，允許以 **$SpO_2/FiO_2 \le 315$** 代替傳統動脈血氣的 $PaO_2/FiO_2 \le 300$。
        3. **影像學標準升級**：
           - 保留雙側肺部浸潤影 (bilateral opacities) 診斷。
           - **新增肺部超音波 (Lung Ultrasound)** 作為可接受的影像診斷工具，極佳適用於不便搬動病患或資源有限環境。
        4. **資源有限地區的簡化條件**：
           - 在這類環境下，取消對 PEEP、特定氧氣流速或特定呼吸器設備的強制要求，降低床邊診斷門檻。
        """, unsafe_allow_html=True)
        
        # Try to show image
        import os
        if os.path.exists("Ards define.jpg"):
            st.image("Ards define.jpg", caption="ARDS 定義演進歷程 (Evolution of ARDS Definition)", use_container_width=True)
        elif os.path.exists("/workspace/Ards define.jpg"):
            st.image("/workspace/Ards define.jpg", caption="ARDS 定義演進歷程 (Evolution of ARDS Definition)", use_container_width=True)
        else:
            st.info("💡 **提示**：若要在網頁端顯示 ARDS 定義演進圖表，可將您下載的 `Ards define.jpg` 放置於本程式相同的路徑下。")
            
        # Severity Table
        st.markdown(r"""
        ### 📊 ARDS 嚴重度分級與全球死亡率：
        | 嚴重度分級 | $PaO_2/FiO_2$ (mmHg) | $SpO_2/FiO_2$ | 全球預估住院死亡率 |
        | :--- | :---: | :---: | :---: |
        | **🟢 輕度 (Mild ARDS)** | 200 ~ 300 | 235 ~ 315 | **35%** |
        | **🟡 中度 (Moderate ARDS)** | 100 ~ 200 | 148 ~ 235 | **40%** |
        | **🔴 重度 (Severe ARDS)** | $\le$ 100 | $\le$ 148 | **46%** |
        """, unsafe_allow_html=True)
        
        st.markdown("---")
        
        # Section B: Personalized Care Strategy
        st.markdown(r"""
        ### 🎯 ARDS 臨床個人化生理照護策略
        現代 ARDS 照護強調整體與**「個人化生理導向 (Individualized Care)」**的結合，擺脫傳統一體適用 (one-size-fits-all) 的限制：
        
        1. **早期病因識別**：
           - 積極清除與治療肺炎、敗血症、誤吸等原發病因。
        2. **臨床生理個人化分類 (計算器對接)**：
           - 📐 **影像形態學**：**局部型 (Focal)** vs. **非局部/彌漫型 (Non-focal)**。局部型應避免高 PEEP，防止正常肺泡過度膨脹。
           - 🧬 **生物表型**：**低發炎型 (Type 1)** vs. **高發炎型 (Type 2)**。Type 2 高發炎患者對 Simvastatin 及類固醇等精準治療反應佳 (可於 Tab 1 計算)。
           - 🔄 **可復張性**：**高復張型 (R/I $\ge$ 0.5)** vs. **低復張型 (R/I < 0.5)**。用以決定病患適合中高 PEEP 還是限制為中低 PEEP (可於 Tab 3 計算)。
        3. **階梯式呼吸支持策略**：
           - **HFNO**：輕度 ARDS 優先，舒適度高，利於咳痰與溝通。
           - **NIV 非侵入正壓**：輕中度無休克患者，優先推薦**頭盔式 (Helmet) 介面**。需於 1-2 小時內評慢評估成效，避免因 NIV 失敗而延誤插管。
           - **IMV 侵入式通氣**：中重度患者首選。插管後先給予安全預設 ($V_t$ 6 mL/kg PBW, 安全 PEEP)，1小時內抽血量測力學並展開個人化調校。
        4. **右心室 (RV) 保護與監測**：
           - 20-25% 患者會併發右心功能不全 (ACP)。高危患者應降低超音波評估門檻，嚴格限制 $P_{plat} < 26\text{–}28 \text{ cmH}_2\text{O}$。
        5. **保守水分管理 (FACCT 策略)**：
           - 休克緩解後採取保守水分限制或利尿策略，追求中性累積水分平衡，能有效減輕肺水腫，顯著縮短 ICU 住院與呼吸器使用天數。
        """, unsafe_allow_html=True)

        st.markdown("---")
        st.markdown(r"""
        ### 📚 核心參考文獻 (Key Reference Papers)
        1. **Wongtirawit N, Menga LS, Brito R, Docci M, Plens GM, Alcala G, Cantan B, Bosma KJ, Ko M, Roman-Sarita G, Purcell T, Cominesi DR, Greendyk RA, Brochard L.**  
           *ARDS management beyond the guidelines: a practical physiology-based approach to individualized care.*  
           **Intensive Care Medicine** (2026). DOI: [10.1007/s00134-026-08563-7](https://doi.org/10.1007/s00134-026-08563-7)
        2. **Morris IS, Amato M, Kassis EB, Bellani G, Calfee CS, Heunks L, Hodgson C, Nair P, Serpa Neto A, Sahetya S, Summers C, Telias I, Yoshida T, Slutsky AS, Ferguson ND.**  
           *The medical management of acute respiratory distress syndrome.*  
           **Intensive Care Medicine** (2026) 52:104–117. DOI: [10.1007/s00134-025-08251-y](https://doi.org/10.1007/s00134-025-08251-y)
        """, unsafe_allow_html=True)


    # ---------------------------------------------------------
    # TAB 1: Biological Phenotypes (Sinha 2020) - MOVED TO TAB 1
    # ---------------------------------------------------------
    with tab1:
        st.markdown('<div class="tab-content"></div>', unsafe_allow_html=True)
        st.subheader("🧬 ARDS 臨床生物發炎亞型預測 (Sinha 2020)")
        st.markdown("""
        根據 Sinha 2020 年發表於《Lancet Respiratory Medicine》的標誌性研究，ARDS 可分為**「低發炎型 (Type 1)」**與**「高發炎型 (Type 2)」**。
        本工具已內建論文中發布的**「黃金主模型 (Primary Models)」**與**「臨床輔助模型 (Ancillary Models)」**。
        """)
        
        model_type = st.radio("選擇預測指標模型", [
            "🏆 Sinha 2020 論文主模型：IL-8 & Protein C 4變數模型 (最推薦)",
            "🏆 Sinha 2020 論文主模型：IL-8 & Protein C 3變數模型 (免收縮劑)",
            "🧬 Sinha 2020 輔助模型：IL-6 & Protein C 4變數模型 (Model 5)",
            "🧬 Sinha 2020 輔助模型：IL-8 & sTNFR-1 4變數模型 (Model 2)"
        ], key="model_selection_t1")
        
        st.markdown("---")
        
        # Form inputs
        col_bio1, col_bio2 = st.columns(2)
        prob = 0.0
        
        if "論文主模型：IL-8 & Protein C 4變數" in model_type:
            with col_bio1:
                il8 = st.number_input("白介素-8 IL-8 (pg/mL)", min_value=1.0, max_value=20000.0, value=30.0, step=5.0, help="低發炎型中位數為 23，高發炎型為 133", key="il8_prim4_t1")
                protein_c = st.number_input("蛋白質 C 活性 Protein C (% control)", min_value=1.0, max_value=250.0, value=80.0, step=5.0, help="低發炎型平均為 96.0，高發炎型為 53.5", key="pc_prim4_t1")
            with col_bio2:
                bicarbonate_bio = st.number_input("碳酸氫根 Bicarbonate (mmol/L)", min_value=2.0, max_value=50.0, value=22.0, step=1.0, help="低發炎型平均為 23.1，高發炎型為 17.3", key="bicarb_prim4_t1")
                vasopressor_bio = st.selectbox("是否使用血管收縮劑 (Vasopressor)", ["否 (No)", "是 (Yes)"], help="高發炎型患者 66% 需使用收縮劑，低發炎型僅 18%", key="vaso_prim4_t1")
                
            vp_bool = True if "是" in vasopressor_bio else False
            prob = calculate_primary_il8_4var_prob(il8, protein_c, bicarbonate_bio, vp_bool)
            
            # Show mathematical formula explicitly with prepended 'r'
            st.markdown("##### 📐 模型計算公式 (Lancet RM 2020 - Primary 4-var):")
            st.markdown(r"""
            $$
            \ln\left(\frac{P}{1-P}\right) = 4.8127 + 1.5424 \times \ln(IL8 + 1) - 0.2502 \times Bicarbonate - 1.8778 \times \ln(Protein C + 1) + 2.1026 \times Vasopressor
            $$
            """, unsafe_allow_html=True)
            
        elif "論文主模型：IL-8 & Protein C 3變數" in model_type:
            with col_bio1:
                il8 = st.number_input("白介素-8 IL-8 (pg/mL)", min_value=1.0, max_value=20000.0, value=30.0, step=5.0, help="低發炎型中位數為 23，高發炎型為 133", key="il8_prim3_t1")
                protein_c = st.number_input("蛋白質 C 活性 Protein C (% control)", min_value=1.0, max_value=250.0, value=80.0, step=5.0, help="低發炎型平均為 96.0，高發炎型為 53.5", key="pc_prim3_t1")
            with col_bio2:
                bicarbonate_bio = st.number_input("碳酸氫根 Bicarbonate (mmol/L)", min_value=2.0, max_value=50.0, value=22.0, step=1.0, help="低發炎型平均為 23.1，高發炎型為 17.3", key="bicarb_prim3_t1")
                st.markdown("<br><p style='font-size:12px; color:#7f8c8d; font-style:italic;'>💡 3變數模型已排除收縮劑影響，更具備跨機構決策穩定度。</p>", unsafe_allow_html=True)
                
            prob = calculate_primary_il8_3var_prob(il8, protein_c, bicarbonate_bio)
            
            # Show mathematical formula explicitly
            st.markdown("##### 📐 模型計算公式 (Lancet RM 2020 - Primary 3-var):")
            st.markdown(r"""
            $$
            \ln\left(\frac{P}{1-P}\right) = 6.1241 + 1.4226 \times \ln(IL8 + 1) - 0.2596 \times Bicarbonate - 1.8330 \times \ln(Protein C + 1)
            $$
            """, unsafe_allow_html=True)
            
        elif "IL-6 & Protein C 4變數" in model_type:
            with col_bio1:
                il6 = st.number_input("白介素-6 IL-6 (pg/mL)", min_value=1.0, max_value=100000.0, value=150.0, step=10.0, help="低發炎型中位數為 116，高發炎型為 933", key="il6_anc4_t1")
                protein_c = st.number_input("蛋白質 C 活性 Protein C (% control)", min_value=1.0, max_value=250.0, value=80.0, step=5.0, help="低發炎型平均為 96.0，高發炎型為 53.5", key="pc_anc4_t1")
            with col_bio2:
                bicarbonate_bio = st.number_input("碳酸氫根 Bicarbonate (mmol/L)", min_value=2.0, max_value=50.0, value=22.0, step=1.0, help="低發炎型平均為 23.1，高發炎型為 17.3", key="bicarb_anc4_t1")
                vasopressor_bio = st.selectbox("是否使用血管收縮劑 (Vasopressor)", ["否 (No)", "是 (Yes)"], key="vaso_anc4_t1")
                
            vp_bool = True if "是" in vasopressor_bio else False
            prob = calculate_ancillary_il6_4var_prob(il6, protein_c, bicarbonate_bio, vp_bool)
            
            # Show mathematical formula explicitly
            st.markdown("##### 📐 模型計算公式 (Lancet RM 2020 - Ancillary Model 5):")
            st.markdown(r"""
            $$
            \ln\left(\frac{P}{1-P}\right) = 4.0323 + 0.9191 \times \ln(IL6 + 1) - 0.2581 \times Bicarbonate - 1.3805 \times \ln(Protein C + 1) + 1.7412 \times Vasopressor
            $$
            """, unsafe_allow_html=True)
            
        else: # IL-8 & sTNFR-1
            with col_bio1:
                il8 = st.number_input("白介素-8 IL-8 (pg/mL)", min_value=1.0, max_value=20000.0, value=30.0, step=5.0, help="低發炎型中位數為 23，高發炎型為 133", key="il8_anc2_t1")
                stnfr1 = st.number_input("可溶性腫瘤壞死因子受體-1 sTNFR-1 (pg/mL)", min_value=100.0, max_value=100000.0, value=3500.0, step=100.0, help="低發炎型中位數為 3225，高發炎型為 7452", key="stnfr_anc2_t1")
            with col_bio2:
                bicarbonate_bio = st.number_input("碳酸氫根 Bicarbonate (mmol/L)", min_value=2.0, max_value=50.0, value=22.0, step=1.0, help="低發炎型平均為 23.1，高發炎型為 17.3", key="bicarb_anc2_t1")
                vasopressor_bio = st.selectbox("是否使用血管收縮劑 (Vasopressor)", ["否 (No)", "是 (Yes)"], key="vaso_anc2_t1")
                
            vp_bool = True if "是" in vasopressor_bio else False
            prob = calculate_ancillary_il8_sTNFR1_4var_prob(il8, stnfr1, bicarbonate_bio, vp_bool)
            
            # Show mathematical formula explicitly
            st.markdown("##### 📐 模型計算公式 (Lancet RM 2020 - Ancillary Model 2):")
            st.markdown(r"""
            $$
            \ln\left(\frac{P}{1-P}\right) = -13.1351 + 1.3947 \times \ln(IL8 + 1) - 0.2145 \times Bicarbonate + 1.1818 \times \ln(sTNFR1 + 1) + 2.1398 \times Vasopressor
            $$
            """, unsafe_allow_html=True)
            
        # Display Phenotype Prediction Results
        st.subheader("📊 亞型預測結果")
        prob_pct = prob * 100.0
        
        st.markdown(f"""
        <div class="metric-container" style="background-color: #f1f2f6;">
            <div class="metric-label">高發炎亞型 (Type 2 / Hyper-inflammatory) 預估機率</div>
            <div class="metric-value" style="color: {'#e74c3c' if prob >= 0.5 else '#2ebd59'};">{prob_pct:.1f} %</div>
            <div style="font-size: 11px; color:#7f8c8d;">臨床判定切點為 50.0%</div>
        </div>
        """, unsafe_allow_html=True)
        
        if prob >= 0.5:
            st.markdown(f"""
            <div class="warning-box" style="background-color: #e74c3c;">
                <h3 style="margin-top:0; color:white; font-size:18px;">🧬 預測亞型：Type 2 高發炎型 (Hyper-inflammatory)</h3>
                <p style="margin-bottom:8px; font-size:14px; font-weight:500;">此患者展現出強烈的全身性發炎反應，預後風險極高 (死亡率與器官衰竭天數顯著高於 Type 1)。</p>
                <hr style="border-top: 1px solid rgba(255,255,255,0.3); margin: 10px 0;">
                <p style="margin-bottom:0; font-size:12px; line-height:1.5;">
                    <strong>💡 「預測性富集」精準醫學決策價值：</strong><br>
                    1. <strong>辛伐他汀 (Simvastatin) 治療反應極佳</strong>：HARP-2 試驗事後分析證使，雖然全體分析為陰性，但 <strong>Type 2 高發炎型患者接受 Simvastatin 治療能顯著提升 28 天存活率</strong>。<br>
                    2. <strong>水分管理可能耐受寬鬆</strong>：FACCT 試驗事後分析顯示，相較於 Type 1 對保守水分有良好反應，Type 2 患者在特定情況下對寬鬆補液可能呈現更好的生存趨勢。<br>
                    3. <strong>高度右心衰竭風險</strong>：應極度嚴防急性肺心症 (ACP)。
                </p>
            </div>
            """, unsafe_allow_html=True)
        else:
            st.markdown(f"""
            <div class="warning-box" style="background-color: #2ebd59;">
                <h3 style="margin-top:0; color:white; font-size:18px;">🧬 預測亞型：Type 1 低發炎型 (Hypoinflammatory)</h3>
                <p style="margin-bottom:8px; font-size:14px; font-weight:500;">全身性發炎反應較為輕微，預後顯著優於高發炎型 (約佔全體 ARDS 患者的 70%)。</p>
                <hr style="border-top: 1px solid rgba(255,255,255,0.3); margin: 10px 0;">
                <p style="margin-bottom:0; font-size:12px; line-height:1.5;">
                    <strong>💡 「預測性富集」精準醫學決策價值：</strong><br>
                    1. <strong>水分限制策略 (FACCT 策略) 反應極佳</strong>：保守限制水分、追求中性累積水分平衡，對此表型能帶來最顯著的減少呼吸器天數與減輕肺水腫效益。<br>
                    2. <strong>不建議常規使用 Simvastatin</strong>：HARP-2 事後分析顯示 Simvastatin 對 Type 1 存活無益，甚至可能有副作用風險。<br>
                    3. <strong>維持穩定保護通氣</strong>：按常規指引維持低潮氣容積、個人化 PEEP，其餘器官衰竭發生機率較低。
                </p>
            </div>
            """, unsafe_allow_html=True)
            
        st.subheader("📚 生物亞型背後的研究背景 (Sinha 2020)")
        st.markdown("""
        *   **什麼是 ARDS 亞型 (Subphenotypes)？**  
            傳統上 ARDS 被視為單一綜合症。但 2014 年 Calfee 團隊利用潛在類別分析 (LCA)，首次在多個隨機對照試驗中鑑定出兩大具有不同發炎和免疫特性的亞型。
        *   **簡化模型 (Parsimonious Models) 的開發意義**：  
            研究最初使用包含 30 餘種臨床 and 分子指標的 LCA 分析。為了實用化，Sinha 2020 年開發了這套「簡化模型」——透過 3-4 個關鍵指標（IL-8, Bicarbonate, Protein C, Vasopressor），可在不損失預測力 (AUC 達 0.94 - 0.96) 的前提下在床邊快速分類，對未來 ARDS 的精準用藥 (Targeted Therapy) 具有奠基地位。
        """)

    # ---------------------------------------------------------
    # TAB 2: Mechanics, Efficiency, PBW Target Grid - MOVED TO TAB 2
    # ---------------------------------------------------------
    with tab2:
        st.markdown('<div class="tab-content"></div>', unsafe_allow_html=True)
        st.subheader("📋 輸入病患基本資料")
        
        input_method = st.radio("PBW (預估體重) 輸入方式", ["由身高性別計算", "直接手動輸入"], key="input_method_t2")
        
        pbw = 60.0
        if input_method == "由身高性別計算":
            col_sex, col_h = st.columns(2)
            with col_sex:
                sex = st.selectbox("病患性別", ["男 (Male)", "女 (Female)"], key="sex_t2")
            with col_h:
                height = st.slider("病患身高 (cm)", 120, 210, 165, key="height_t2")
            sex_str = "male" if "男" in sex else "female"
            pbw = calculate_pbw(height, sex_str)
            st.info(f"計算得預估體重 (PBW): **{pbw:.1f} kg**")
        else:
            pbw = st.number_input("預估體重 PBW (kg)", min_value=10.0, max_value=200.0, value=60.0, step=1.0, key="pbw_manual_t2")
            
        # Display PBW formulas mathematically as requested by user
        st.subheader("📐 PBW 預估體重計算公式 (Predicted Body Weight)")
        st.markdown(r"""
        *   **男 (Male)**: $\text{PBW} = 50.0 + 0.91 \times (\text{身高 cm} - 152.4)$
        *   **女 (Female)**: $\text{PBW} = 45.5 + 0.91 \times (\text{身高 cm} - 152.4)$
        """, unsafe_allow_html=True)

        # Display VR & Vd/Vt formulas mathematically and note on predicted PaCO2
        st.subheader("📐 通氣比例 (VR) 與生理死腔估算核心公式 (Sinha 2019 AJRCCM)")
        st.markdown(r"""
        1. **通氣比例 (Ventilatory Ratio, VR) 公式**:
           $$	ext{VR} = rac{V_E (	ext{L/min}) 	imes PaCO_2 (	ext{mmHg})}{	ext{Predicted } V_E 	imes 	ext{Predicted } PaCO_2} = rac{V_E 	imes PaCO_2}{4.0 	imes 	ext{PBW}}$$
        2. **預估生理死腔比例 ($V_d/V_t$) 估算公式**:
           $$V_d/V_t = 1.0 - rac{0.70}{	ext{VR}} \quad (	ext{若 } 	ext{VR} \le 0.7 	ext{ 則設定底限為 } 0.30)$$

        > 💡 **關於預測 $PaCO_2$ ($37.5$ vs. $40 	ext{ mmHg}$) 臨床數學與生理說明**：  
        > *   **理論推導 ($37.5 	ext{ mmHg}$ / 分母 $3.75 	imes 	ext{PBW}$)**：在 Sinha 2019 原始論文之單位轉換推導中，預測動脈 $PaCO_2$ 採用國際單位制正中間值 $5.0 	ext{ kPa} pprox 37.5 	ext{ mmHg}$。搭配預估分鐘通氣量 $	ext{Predicted } V_E = 100 	ext{ mL/kg/min} 	imes 	ext{PBW} = 0.1 	ext{ L/kg/min} 	imes 	ext{PBW}$，理論分母為 $0.1 	imes 37.5 	imes 	ext{PBW} = \mathbf{3.75 	imes 	ext{PBW}}$。  
        > *   **臨床簡化與應用 ($40 	ext{ mmHg}$ / 分母 $4.0 	imes 	ext{PBW}$)**：為了方便床邊心算與貼合臨床慣用理想正常值 $PaCO_2 = 40 	ext{ mmHg}$，臨床上普遍將分母四捨五入簡化為 **$4.0 	imes 	ext{PBW}$**（即 $0.1 	imes 40 	imes 	ext{PBW}$）。兩者計算出的死腔比例極微小，完全不影響風險分級與床邊決策。本計算器採用臨床通用之 **$4.0 	imes 	ext{PBW}$** 標準。
        """, unsafe_allow_html=True)

        st.subheader("🎛️ 輸入呼吸器與力學參數")
        col_res1, col_res2 = st.columns(2)
        with col_res1:
            vt = st.number_input("潮氣容積 Vt (mL)", min_value=100, max_value=1000, value=360, step=10, key="vt_t2")
            pplat = st.number_input("平台壓 Pplat (cmH2O)", min_value=5, max_value=60, value=25, step=1, key="pplat_t2")
            total_peep = st.number_input("總 PEEP (含 intrinsic PEEP, cmH2O)", min_value=0, max_value=40, value=10, step=1, key="peep_t2")
        with col_res2:
            rr = st.number_input("呼吸速率 RR (bpm)", min_value=5, max_value=50, value=20, step=1, key="rr_t2")
            ve = st.number_input("實際分鐘通氣量 Ve (L/min)", min_value=1.0, max_value=40.0, value=12.0, step=0.5, key="ve_t2")
            paco2 = st.number_input("實際動脈血 PaCO2 (mmHg)", min_value=10.0, max_value=150.0, value=50.0, step=1.0, key="paco2_t2")
            
        st.subheader("⚠️ 氣道開啟壓 (AOP) 偵測與校正 (選填)")
        has_aop = st.checkbox("病患有氣道關閉 (Airway Closure)，需啟用 AOP 校正", key="has_aop_t2")
        aop_val = 0.0
        if has_aop:
            aop_val = st.number_input("測得之氣道開啟壓 AOP (cmH2O)", min_value=1.0, max_value=45.0, value=12.0, step=1.0, key="aop_val_t2")
            
        # Standard Calculations
        vr = calculate_vr(ve, paco2, pbw)
        vd_vt = estimate_vd_vt(vr)
        
        # Calculate corrected driving pressure & compliance
        if has_aop and total_peep < aop_val:
            dp = pplat - aop_val
            st.warning(f"💡 PEEP ({total_peep:.1f}) 低於 AOP ({aop_val:.1f})！已自動使用 AOP 校正驅動壓: ΔP_corrected = Pplat - AOP = **{dp:.1f} cmH2O**")
        else:
            dp = pplat - total_peep
            
        compliance = vt / (pplat - total_peep) if (pplat - total_peep) > 0 else 0.0
        vt_per_pbw = vt / pbw if pbw > 0 else 0.0
        vi = 4 * dp + rr # Ventilation Intensity
        
        # Grid layout for outputs
        st.subheader("📈 床邊生理計算結果")
        col_m1, col_m2, col_m3 = st.columns(3)
        with col_m1:
            st.markdown(f"""
            <div class="metric-container">
                <div class="metric-label">通氣比例 (VR)</div>
                <div class="metric-value">{vr:.2f}</div>
                <div style="font-size: 11px; color:#7f8c8d;">正常人為 1.0</div>
            </div>
            """, unsafe_allow_html=True)
            st.markdown(f"""
            <div class="metric-container">
                <div class="metric-label">通氣強度 (VI)</div>
                <div class="metric-value">{vi:.1f}</div>
                <div style="font-size: 11px; color:#7f8c8d;">公式: 4×ΔP + RR</div>
            </div>
            """, unsafe_allow_html=True)
        with col_m2:
            st.markdown(f"""
            <div class="metric-container">
                <div class="metric-label">預估死腔比例 (Vd/Vt)</div>
                <div class="metric-value">{vd_vt:.2f}</div>
                <div style="font-size: 11px; color:#7f8c8d;">健康上限 0.30</div>
            </div>
            """, unsafe_allow_html=True)
            st.markdown(f"""
            <div class="metric-container">
                <div class="metric-label">靜態順應性 (Crs)</div>
                <div class="metric-value">{compliance:.1f}</div>
                <div style="font-size: 11px; color:#7f8c8d;">mL/cmH2O</div>
            </div>
            """, unsafe_allow_html=True)
        with col_m3:
            st.markdown(f"""
            <div class="metric-container">
                <div class="metric-label">實測驅動壓 (ΔP)</div>
                <div class="metric-value" style="color: {'#e74c3c' if dp >= 15 else '#2ebd59'};">{dp:.1f}</div>
                <div style="font-size: 11px; color:#7f8c8d;">安全切點 < 15</div>
            </div>
            """, unsafe_allow_html=True)
            st.markdown(f"""
            <div class="metric-container">
                <div class="metric-label">實測潮氣/PBW</div>
                <div class="metric-value">{vt_per_pbw:.1f}</div>
                <div style="font-size: 11px; color:#7f8c8d;">mL/kg PBW</div>
            </div>
            """, unsafe_allow_html=True)
            
        # Decsison Box for Tab 2
        risk = get_risk_status(vr, vd_vt)
        st.markdown(f"""
        <div class="warning-box" style="background-color: {risk['color']};">
            <h4 style="margin-top:0; color:white;">🚨 死腔與預後警訊：{risk['title']}</h4>
            <p style="margin-bottom:8px; font-size:13px;">{risk['desc']}</p>
            <p style="margin-bottom:0; font-size:12px; line-height:1.4;"><strong>臨床決策建議：</strong><br>{risk['action'].replace('\\n', '<br>')}</p>
        </div>
        """, unsafe_allow_html=True)
        
        # VR & Vd/Vt Reference ranges card directly visible under the decision block (as requested)
        st.subheader("📊 VR & $V_d/V_t$ 臨床危險區間參考範圍 (Risk Reference Card)")
        st.markdown(r"""
        | 警訊分級 | 通氣比例 (VR) 範圍 | 預估生理死腔比例 ($V_d/V_t$) 範圍 | 生理狀態與 Nuckton 2002 臨床預後意義 |
        | :---: | :---: | :---: | :--- |
        | **🟢 低風險區** | **VR < 1.50** | **$V_d/V_t$ < 0.54** | 通氣效率優良。低於 ARDS 存活患者之平均死腔比值。 |
        | **🟡 中度警戒區** | **1.50 ~ 1.90** | **0.54 ~ 0.63** | 通氣效率受損。介於存活組與死亡組平均值之間，死亡風險隨死腔增加而遞增。 |
        | **🔴 高度危險區** | **1.90 ~ 2.25** | **0.63 ~ 0.69** | 通氣效率嚴重受損。已達或超過 ARDS 死亡組的平均死腔比值 ($V_d/V_t \ge 0.63$)。 |
        | **🛑 極高度危殆** | **VR $\ge$ 2.25** | **$V_d/V_t \ge$ 0.69** | 死腔比例極高。落入 Nuckton 研究中最嚴重的前 20% (Quintile 5) 極高危死亡人群。 |
        """, unsafe_allow_html=True)

        # Target Tidal Volume Grid
        st.subheader("🎯 臨床保護性潮氣容積對照表")
        st.markdown(f"""
        為了落實保護性肺通氣 (Lung-protective Ventilation)，請將病患的潮氣容積鎖定在 **4 至 8 mL/kg PBW**：
        - 🟢 **6 mL/kg PBW (黃金標準)**: **{6 * pbw:.0f} mL**
        - 🟢 **5 mL/kg PBW (嚴重順應性差)**: **{5 * pbw:.0f} mL**
        - 🟡 **4 mL/kg PBW (極致超保護肺通氣)**: **{4 * pbw:.0f} mL** *(此時應使用加熱加濕器，避免使用人工鼻以減少儀器死腔)*
        - 🟡 **8 mL/kg PBW (高順應性放寬上限)**: **{8 * pbw:.0f} mL**
        """, unsafe_allow_html=True)


    # ---------------------------------------------------------
    # TAB 3: Lung Recruitability (R/I Ratio) - MOVED TO TAB 3
    # ---------------------------------------------------------
    with tab3:
        st.markdown('<div class="tab-content"></div>', unsafe_allow_html=True)
        st.subheader("🔄 肺可復張性評估 (Recruitment-to-Inflation Ratio)")
        st.markdown("""
        根據 Wongtirawit 2026 (ICM) 綜述，我們不應對所有 ARDS 病患無差別地調高 PEEP。\n
        使用 **R/I Ratio 呼氣單步降壓法**，可以幫助醫師在床邊 10 秒內區分出病患屬於 **「高復張性」** 還是 **「低復張性」** 亞型，避免造成嚴重肺部過度充氣與右心後負荷增加。\n
        """)
        
        # Display mathematical formulas for R/I Ratio as requested by user
        st.subheader("📐 R/I Ratio 核心計算生理公式")
        st.markdown(r"""
        1. **低 PEEP 靜態順應性**: 
           $$C_{\text{low}} = \frac{Vt_{\text{low}}}{P_{\text{plat,low}} - PEEP_{\text{low}}}$$
        2. **單步降壓之總呼氣末肺容積變化量**: 
           $$\Delta EELV = V_{\text{exp,total}} - Vt_{\text{high}}$$
        3. **預期未復張之充氣肺泡容積**: 
           $$V_{\text{predicted}} = C_{\text{low}} \times (PEEP_{\text{high}} - PEEP_{\text{low}})$$
        4. **高 PEEP 真正復張(救回)的肺泡容積**: 
           $$V_{\text{recruited}} = \Delta EELV - V_{\text{predicted}}$$
        5. **最終 R/I Ratio 比值**: 
           $$R/I \text{ Ratio} = \frac{V_{\text{recruited}}}{(PEEP_{\text{high}} - PEEP_{\text{low}}) \times C_{\text{low}}}$$
        """, unsafe_allow_html=True)

        st.subheader("📋 R/I Ratio 床邊詳細操作步驟 (Step-by-Step Bedside Guide)")
        st.markdown(r"""
        <div style="background-color: #f0f7fc; padding: 18px; border-radius: 8px; border-left: 5px solid #2980b9; margin-bottom: 20px;">
            <p style="margin-top:0; font-size:15px; font-weight:600; color:#2980b9;">🛠️ 呼氣單步降壓法 (One-Breath PEEP Reduction) 執行步驟：</p>
            <ol style="margin-bottom:8px; font-size:13px; line-height:1.6; color:#2c3e50; padding-left:20px;">
                <li><strong>高壓穩定階段 (PEEPhigh Stabilization)</strong>：將病患呼吸器之 PEEP 設為高水平 (通常設定為 <b>15 cmH₂O</b>)，並穩定通氣 <b>10 分鐘</b>，使肺泡充分復張並達到穩態。</li>
                <li><strong>防氣體陷縮 & 記錄基準潮氣量</strong>：短暫將呼吸速率 (RR) 調降至 <b>6–8 bpm</b> (進行 1-2 次呼吸)，以完全排除氣道捕獲氣體 (Air trapping / Auto-PEEP)。記錄此時的呼出潮氣容積 <b>Vt_high</b> ($VTe_{\text{high}}$，如 450 mL)。</li>
                <li><strong>單步突發降壓 (One-Breath Drop)</strong>：在<b>單次吐氣開始時</b>，將 PEEP 瞬間調降至低水平 (通常設定為 <b>5 cmH₂O</b>)。此時病患會因為 sudden drop 而吐出一口極大的氣體。</li>
                <li><strong>記錄總呼出容積 (V_exp_drop)</strong>：在降壓那一瞬間的單個呼吸週期中，記錄呼吸器面板上測得的<b>總呼出容積</b> (V_exp_total / V_exp_drop，如 1100 mL)。</li>
                <li><strong>低壓生理力學量測 (PEEPlow Mechanics)</strong>：讓患者在低 PEEP 水平穩定呼吸 2-3 次。短暫按下「吸氣阻斷鍵 (Inspiratory Hold) 0.2-0.3 秒」測量低 PEEP 下的平台壓 <b>Pplat_low</b> 與呼出潮氣 <b>Vt_low</b>。</li>
                <li><strong>輸入參數計算</strong>：將上述測得之 6 個數值代入下方對應之輸入欄位中，一秒取得您的 R/I Ratio 報告！</li>
            </ol>
            <hr style="border-top: 1px solid rgba(41, 128, 185, 0.2); margin: 10px 0;">
            <p style="margin-bottom:0; font-size:12px; line-height:1.5; color:#2c3e50;">
                ⚠️ <b>AOP 關鍵提醒 (Airway Opening Pressure)</b>：<br>
                在執行 R/I 測試前，強烈建議先進行低流速吸氣量測 AOP。如果病患的 <b>AOP 高於 PEEPlow</b> (例如 AOP 為 8 cmH₂O，而您設定之 PEEPlow 為 5 cmH₂O)，<b>您必須使用 AOP 作為 PEEPlow</b> 代替 5 cmH₂O！否則，在降壓時氣道會提前完全關閉，導致順應性與 R/I Ratio 算得嚴重失真。
            </p>
        </div>
        """, unsafe_allow_html=True)

        col_ri1, col_ri2 = st.columns(2)
        with col_ri1:
            st.markdown("**1. PEEP 降壓前參數設定**")
            peep_high = st.number_input("高 PEEP 水平 PEEPhigh (cmH2O)", min_value=10, max_value=30, value=15, step=1, key="peep_high_t2")
            vt_high = st.number_input("PEEPhigh 下的呼出潮氣容積 VTe (mL)", min_value=100, max_value=1000, value=450, step=10, key="vte_high_t2")
            v_exp_drop = st.number_input("一步降壓至 PEEPlow 時，排出之「總呼出容積」(mL)", min_value=200, max_value=2500, value=1100, step=50, key="v_exp_drop_t2", help="指呼吸器在 PEEP 瞬間降壓後那一呼吸週期的總呼出容積")
        with col_ri2:
            st.markdown("**2. PEEP 降壓後生理力學**")
            peep_low = st.number_input("低 PEEP 水平 PEEPlow (cmH2O)", min_value=0, max_value=15, value=5, step=1, key="peep_low_t2")
            vt_low = st.number_input("PEEPlow 下的潮氣容積 Vt (mL)", min_value=100, max_value=1000, value=400, step=10, key="vt_low_t2")
            pplat_low = st.number_input("PEEPlow 下的平台壓 Pplat (cmH2O)", min_value=5, max_value=50, value=15, step=1, key="pplat_low_t2")
            
        ri_ratio, c_low, delta_eelv, v_recruited = calculate_ri_ratio(peep_high, peep_low, vt_high, v_exp_drop, pplat_low, vt_low)
        
        st.subheader("📊 R/I Ratio 評估結果")
        col_ri_res1, col_ri_res2, col_ri_res3 = st.columns(3)
        with col_ri_res1:
            st.markdown(f"""
            <div class="metric-container">
                <div class="metric-label">R/I Ratio</div>
                <div class="metric-value" style="color: {'#2ebd59' if ri_ratio >= 0.5 else '#e74c3c'};">{ri_ratio:.2f}</div>
                <div style="font-size: 11px; color:#7f8c8d;">指標切點 0.5</div>
            </div>
            """, unsafe_allow_html=True)
        with col_ri_res2:
            st.markdown(f"""
            <div class="metric-container">
                <div class="metric-label">復張肺容積 (V_recruited)</div>
                <div class="metric-value">{v_recruited:.0f} mL</div>
                <div style="font-size: 11px; color:#7f8c8d;">高 PEEP 救回的肺泡量</div>
            </div>
            """, unsafe_allow_html=True)
        with col_ri_res3:
            st.markdown(f"""
            <div class="metric-container">
                <div class="metric-label">EELV 改变量 (ΔEELV)</div>
                <div class="metric-value">{delta_eelv:.0f} mL</div>
                <div style="font-size: 11px; color:#7f8c8d;">呼氣末總肺容積變化</div>
            </div>
            """, unsafe_allow_html=True)
            
        if ri_ratio >= 0.5:
            st.markdown(f"""
            <div class="warning-box" style="background-color: #2ebd59;">
                <h4 style="margin-top:0; color:white;">🟢 高可復張性 (High Recruiter - R/I ≥ 0.5)</h4>
                <p style="margin-bottom:8px; font-size:13px;">該患者的肺部在調高 PEEP 時，能成功救回大量的塌陷肺泡 (V_recruited: {v_recruited:.0f} mL)，增加可參與通氣的工作面積，並能解除缺氧性肺血管收縮。</p>
                <hr style="border-top: 1px solid rgba(255,255,255,0.3); margin: 8px 0;">
                <p style="margin-bottom:0; font-size:12px; line-height:1.4;">
                    <strong>臨床決策導航：</strong><br>
                    1. <strong>適用中高 PEEP 策略</strong> (如 12-20 cmH2O)，能有效防範 Atelectrauma。<br>
                    2. <strong>建議滴定 PEEP</strong>：可進一步結合食道壓量測 (TPP 靶值 0 ± 2 cmH2O) 或電阻抗斷層成像 (EIT) 尋找兼顧肺泡塌陷與過度膨脹的最佳 PEEP。
                </p>
            </div>
            """, unsafe_allow_html=True)
        else:
            st.markdown(f"""
            <div class="warning-box" style="background-color: #e74c3c;">
                <h4 style="margin-top:0; color:white;">🚨 低可復張性 (Low Recruiter - R/I < 0.5)</h4>
                <p style="margin-bottom:8px; font-size:13px;">該患者調高 PEEP 幾乎無法打開新肺泡。此時強行調高 PEEP，多餘的壓力「只會吹破健康肺泡」引發過度膨脹，並壓迫微血管、暴增肺血管阻力 (PVR) 引發右心衰竭。</p>
                <hr style="border-top: 1px solid rgba(255,255,255,0.3); margin: 8px 0;">
                <p style="margin-bottom:0; font-size:12px; line-height:1.4;">
                    <strong>臨床決策導航：</strong><br>
                    1. <strong>限制為中低 PEEP</strong> (設定在 8 - 12 cmH2O 之間即可)。<br>
                    2. <strong>絕對避免激進的肺泡復張手法 (Recruitment Maneuvers, RM)</strong>，避免氣壓傷與循環休克。<br>
                    3. 應優先評估是否存在局部性病灶 (Focal ARDS)，並嚴防右心衰竭 (ACP)。
                </p>
            </div>
            """, unsafe_allow_html=True)


    # ---------------------------------------------------------
    # TAB 4: Spontaneous Breathing Drive & Effort - MOVED TO TAB 4
    # ---------------------------------------------------------
    with tab4:
        st.markdown('<div class="tab-content"></div>', unsafe_allow_html=True)
        st.subheader("🧠 自主呼吸 Efforts 與驅力監測")
        st.markdown("""
        當 ARDS 患者從「控制通氣」過渡到「輔助/自主通氣 (Assisted Ventilation)」時，病患強烈的自主呼吸努力會對肺部造成隱形、不可測的物理拉扯。\n
        這會導致 **病患自殘性肺損傷 (PSILI, Patient Self-Inflicted Lung Injury)** 並且惡化肺水腫。\n
        本工具引進了 Wongtirawit 2026 (ICM) 推薦的床邊無創遮斷指標，防範 PSILI。\n
        """)
        
        # Present safe vs unsafe clinical ranges for assisted ventilation (as requested by user)
        st.subheader("📊 自主輔助通氣下床邊生理指標安全區間 (Spontaneous Effort Safety Ranges)")
        st.markdown(r"""
        | 生理監測指標 | 🟢 正常/安全區 (Safe Range) | 🚨 高危/病損區 (Unsafe / PSILI Risk) | 臨床床邊生理學意義 |
        | :--- | :---: | :---: | :--- |
        | **P0.1 (氣道百毫秒遮斷壓)** | **< 3.5 cmH₂O** | **$\ge$ 3.5 ~ 4.0 cmH₂O** | 反映大腦中樞呼吸驅力（Respiratory Drive）。過高代表患者通氣竭力。 |
        | **$\Delta P_{occ}$ (遮斷努力壓差)** | **$\ge$ -20 cmH₂O** | **< -20 cmH₂O** (如 -25) | 用於估算動態跨肺驅動壓。更負的數值提示橫膈吸氣努力拉扯過度。 |
        | **PMI (肌肉壓力指數)** | **$\le$ 3.0 cmH₂O** | **> 3.0 cmH₂O** | 計算: $P_{\text{plat}} - P_{\text{peak}}$。大於 3 代表病患呼吸肌做功過大，承受高剪力。 |
        """, unsafe_allow_html=True)

        col_drv1, col_res_drv = st.columns(2)
        with col_drv1:
            p01 = st.number_input("P0.1 氣道遮斷壓 (cmH2O)", min_value=0.0, max_value=15.0, value=2.0, step=0.5, help="吸氣開始 100 毫秒時的壓力降，反映大腦呼吸驅力。正常人或合適通氣支持下約為 1.5 - 3.5")
            dp_occ = st.number_input("ΔPocc 呼氣末遮斷努力壓差 (cmH2O)", min_value=-50.0, max_value=0.0, value=-12.0, step=1.0, help="進行呼氣末遮斷(end-expiratory hold)時，病患自主吸氣造成的最深壓力降，用以估算動態跨肺驅動壓")
            
            p_peak_as = st.number_input("自主輔助通氣下 Peak Pressure (cmH2O)", min_value=5.0, max_value=50.0, value=18.0, step=1.0, key="pp_t3")
            p_peak_as = p_peak_as  # dummy assignment
            p_plat_as = st.number_input("自主輔助通氣下 Plateau Pressure (cmH2O)", min_value=5.0, max_value=50.0, value=22.0, step=1.0, key="pplat_t3", help="在 PSV 下短暫按下吸氣遮斷(Inspiratory hold)測得。若因肌肉放鬆導致壓力回彈，則 Pplat 會大於 Ppeak")
            
        pmi = p_plat_as - p_peak_as
        
        with col_res_drv:
            st.markdown("**🧠 呼吸驅力與努力度評估結果**")
            
            # P0.1 evaluation
            p01_status = "🟢 正常/安全呼吸驅力" if p01 < 3.5 else "🚨 驅力過高 (High Drive)"
            p01_color = "#2ebd59" if p01 < 3.5 else "#e74c3c"
            
            # ΔPocc evaluation
            dp_occ_status = "🟢 正常/安全呼吸努力" if dp_occ >= -20.0 else "🚨 努力過度 (Excessive Effort)"
            dp_occ_color = "#2ebd59" if dp_occ >= -20.0 else "#e74c3c"
            
            st.markdown(f"""
            <div class="metric-container" style="background-color: #f8f9fa;">
                <div class="metric-label">P0.1 呼吸驅力</div>
                <div class="metric-value" style="color: {p01_color};">{p01:.1f} cmH2O</div>
                <div style="font-size: 11px; color:#7f8c8d;">安全切點: &lt; 3.5 cmH2O | 狀態: <strong>{p01_status}</strong></div>
            </div>
            """, unsafe_allow_html=True)
            
            st.markdown(f"""
            <div class="metric-container" style="background-color: #f8f9fa;">
                <div class="metric-label">ΔPocc 遮斷努力值</div>
                <div class="metric-value" style="color: {dp_occ_color};">{dp_occ:.1f} cmH2O</div>
                <div style="font-size: 11px; color:#7f8c8d;">安全切點: &ge; -20 cmH2O | 狀態: <strong>{dp_occ_status}</strong></div>
            </div>
            """, unsafe_allow_html=True)
            
            st.markdown(f"""
            <div class="metric-container" style="background-color: #f8f9fa;">
                <div class="metric-label">PMI (Pressure Muscle Index)</div>
                <div class="metric-value">{pmi:.1f} cmH2O</div>
                <div style="font-size: 11px; color:#7f8c8d;">PMI &gt; 3 代表病患自主吸氣做功過大</div>
            </div>
            """, unsafe_allow_html=True)
            
        # Clinical Interventions Recommendation Box for Spontaneous Effort
        if p01 >= 3.5 or dp_occ < -20.0 or pmi > 3.0:
            st.markdown(f"""
            <div class="warning-box" style="background-color: #e74c3c;">
                <h4 style="margin-top:0; color:white;">🚨 高危病患自殘性肺損傷 (PSILI) 預警</h4>
                <p style="margin-bottom:8px; font-size:13px;">病患此時展現出極度強烈的吸氣努力與驅力！這在微觀下會引發 occult pendelluft (肺內氣體異常震盪與重新分配)，加重區域剪力拉扯，引起嚴重的肺水腫。</p>
                <hr style="border-top: 1px solid rgba(255,255,255,0.3); margin: 8px 0;">
                <p style="margin-bottom:0; font-size:12px; line-height:1.4;">
                    <strong>臨床床邊生理調校建議：</strong><br>
                    1. <strong>優化呼吸器設定 (而非直接加深鎮靜)</strong>：調整吸氣流速 (VCV 下增加流速解決 Flow starvation)；調整 PSV 的壓力支持水平。<br>
                    2. <strong>調高 PEEP</strong>：有時肺部 derecruitment 會異常放大吸氣拉扯，調高 PEEP 改善順應性後，呼吸努力反而會下降。<br>
                    3. <strong>控制可逆病因</strong>：積極排除發燒、疼痛、嚴重代謝性酸中毒與二氧化碳滯留等呼吸驅力刺激因子。<br>
                    4. <strong>考慮神經肌肉阻斷劑 (NMBA)</strong>：若經上述調整自主 efforts 仍大於安全上限，應果斷早期使用肌鬆劑癱瘓膈肌，徹底保護肺臟。
                </p>
            </div>
            """, unsafe_allow_html=True)
        else:
            st.markdown(f"""
            <div class="warning-box" style="background-color: #2ebd59;">
                <h4 style="margin-top:0; color:white;">🟢 自主呼吸努力處於安全肺保護區</h4>
                <p style="margin-bottom:8px; font-size:13px;">該病患的 P0.1 與努力做功皆在可控範圍內，此時自主呼吸對橫膈肌具有良好的鍛鍊作用，有助於加速脫離呼吸器。</p>
                <hr style="border-top: 1px solid rgba(255,255,255,0.3); margin: 8px 0;">
                <p style="margin-bottom:0; font-size:12px; line-height:1.4;">
                    <strong>臨床決策：</strong><br>
                    - 維持輕度/每日喚醒鎮靜策略 (Light Sedation)。<br>
                    - 若病患原發病因已好轉，可積極考慮過渡到 proportional 模式 (如 PAV+ / NAVA) 以獲得更完美的自主呼吸人機同步。
                </p>
            </div>
            """, unsafe_allow_html=True)

    st.markdown("---")
    st.caption("聲明：本工具僅供臨床醫學學術討論與生理機制模擬使用，實際呼吸器設定與病人處置應由專科醫師依病患臨床即時狀態做出決定。")

# ---------------------------------------------------------
# CLI Command Line Mode
# ---------------------------------------------------------
def run_cli():
    print("="*65)
    print("      ARDS BEDSIDE PHYSIOLOGICAL & BIOLOGICAL PHENOTYPE MASTER CALCULATOR")
    print("                臨床生理、死腔、可復張性與生物表型計算器 v9")
    print("="*65)
    
    print("\n請選擇您要執行的功能：")
    print(" [1] 🧬 預測 ARDS 生物臨床發炎亞型 (Sinha 2020)")
    print(" [2] 🫁 計算 Ventilatory Ratio (VR) 與生理死腔 (Nuckton 2002 / Sinha 2019)")
    print(" [3] 🔄 計算 Recruitment-to-Inflation (R/I) Ratio (Wongtirawit 2026)")
    print(" [4] 🧠 評估自主呼吸努力與驅力 (P0.1, ΔPocc, PMI)")
    
    choice = input("\n請輸入 1, 2, 3 或 4: ").strip()
    
    if choice == '1':
        print("\n" + "="*50)
        print("          ARDS 生物臨床發炎亞型預測 (Sinha 2020)")
        print("="*50)
        print("請選擇模型類型：")
        print(" [1] 🏆 Sinha 2020 論文主模型：IL-8 & Protein C 4變數模型 (最推薦)")
        print(" [2] 🏆 Sinha 2020 論文主模型：IL-8 & Protein C 3變數模型")
        print(" [3] 🧬 Sinha 2020 輔助模型：IL-6 & Protein C 4變數模型")
        print(" [4] 🧬 Sinha 2020 輔助模型：IL-8 & sTNFR-1 4變數模型")
        m_choice = input("請輸入 1, 2, 3 或 4: ").strip()
        try:
            bicarbonate = float(input("\n請輸入 碳酸氫根 Bicarbonate (mmol/L) [例如 22]: ").strip())
            if m_choice in ['1', '3', '4'] :
                vaso_input = input("是否使用血管收縮劑? (y/n): ").strip().lower()
                vp_bool = True if vaso_input in ['y', 'yes', '是'] else False
            else:
                vp_bool = False
            
            if m_choice == '1':
                il8 = float(input("請輸入 IL-8 (pg/mL) [例如 30]: ").strip())
                protein_c = float(input("請輸入 Protein C (% control) [例如 80]: ").strip())
                prob = calculate_primary_il8_4var_prob(il8, protein_c, bicarbonate, vp_bool)
            elif m_choice == '2':
                il8 = float(input("請輸入 IL-8 (pg/mL) [例如 30]: ").strip())
                protein_c = float(input("請輸入 Protein C (% control) [例如 80]: ").strip())
                prob = calculate_primary_il8_3var_prob(il8, protein_c, bicarbonate)
            elif m_choice == '3':
                il6 = float(input("請輸入 IL-6 (pg/mL) [例如 150]: ").strip())
                protein_c = float(input("請輸入 Protein C (% control) [例如 80]: ").strip())
                prob = calculate_ancillary_il6_4var_prob(il6, protein_c, bicarbonate, vp_bool)
            else:
                il8 = float(input("請輸入 IL-8 (pg/mL) [例如 30]: ").strip())
                stnfr1 = float(input("請輸入 sTNFR-1 (pg/mL) [例如 3500]: ").strip())
                prob = calculate_ancillary_il8_sTNFR1_4var_prob(il8, stnfr1, bicarbonate, vp_bool)
                
            print("-"*50)
            print(f" 高發炎亞型 (Type 2) 機率: {prob*100.0:.1f}%")
            if prob >= 0.5:
                print("🚨 預測亞型：Type 2 高發炎型 (Hyper-inflammatory) - 建議考慮 Simvastatin，防範右心衰竭！")
            else:
                print("🟢 預測亞型：Type 1 低發炎型 (Hypoinflammatory) - 建議採用 FACCT 保守限制輸液。")
            print("="*50 + "\n")
        except ValueError:
            print("\n[錯誤] 輸入數值格式不正確。")
            
    elif choice == '2':
        # Standard VR calculation CLI
        try:
            print("\n[1] 病患預估體重 (PBW) 計算")
            calc_choice = input("是否需要從身高與性別計算 PBW? (y/n, 預設為直接手動輸入): ").strip().lower()
            if calc_choice == 'y':
                sex_input = input("病患性別 (m=男, f=女): ").strip().lower()
                sex = 'male' if sex_input in ['m', 'male', '男'] else 'female'
                height = float(input("病患身高 (cm): ").strip())
                pbw = calculate_pbw(height, sex)
                print(f"--> 計算得出預估體重 (PBW): {pbw:.2f} kg (男: 50 + 0.91*(H-152.4) | 女: 45.5 + 0.91*(H-152.4))")
            else:
                pbw = float(input("請輸入預估體重 PBW (kg): "))
                
            print("\n[2] 呼吸器與血氧參數")
            ve = float(input("實際分鐘通氣量 Ve (L/min): ").strip())
            paco2 = float(input("實際動脈血 PaCO2 (mmHg): ").strip())
            
            vr = calculate_vr(ve, paco2, pbw)
            vd_vt = estimate_vd_vt(vr)
            risk = get_risk_status(vr, vd_vt)
            
            print("\n" + "="*50)
            print("計算結果 (Bedside Results):")
            print(f" - 通氣比例 (VR): {vr:.2f}")
            print(f" - 預估生理死腔比例 (Vd/Vt): {vd_vt:.2f}")
            print("-"*50)
            print(f"🚨 警訊分級: {risk['title']}")
            print(f"📊 生理狀態: {risk['desc']}")
            print(f"👉 臨床決策建議:\n{risk['action']}")
            print("="*50 + "\n")
        except ValueError:
            print("\n[錯誤] 輸入數值格式不正確。")

    elif choice == '3':
        print("\n" + "="*50)
        print("          R/I Ratio 呼氣單步降壓可復張性評估")
        print("="*50)
        try:
            peep_high = float(input("高 PEEP 水平 PEEPhigh (cmH2O) [例如 15]: ").strip())
            vt_high = float(input("PEEPhigh 下的呼出潮氣容積 VTe (mL) [例如 450]: ").strip())
            v_exp_drop = float(input("一步降壓至 PEEPlow 時總呼出容積 (mL) [例如 1100]: ").strip())
            peep_low = float(input("低 PEEP 水平 PEEPlow (cmH2O) [例如 5]: ").strip())
            vt_low = float(input("PEEPlow 下的潮氣容積 Vt (mL) [例如 400]: ").strip())
            pplat_low = float(input("PEEPlow 下的平台壓 Pplat (cmH2O) [例如 15]: ").strip())
            
            ri, c_low, delta_eelv, v_rec = calculate_ri_ratio(peep_high, peep_low, vt_high, v_exp_drop, pplat_low, vt_low)
            print("-" * 50)
            print(" 計算公式對照:")
            print(f"  - Low-PEEP Compliance: C_low = Vt_low / (Pplat_low - PEEP_low) = {c_low:.2f} mL/cmH2O")
            print(f"  - delta_EELV = V_exp_drop - Vt_high = {delta_eelv:.0f} mL")
            print(f"  - V_recruited = delta_EELV - [C_low * (PEEPhigh - PEEPlow)] = {v_rec:.0f} mL")
            print("-" * 50)
            print(f" 計算得出 R/I Ratio: {ri:.2f}")
            print(f" 復張肺容積 (V_recruited): {v_rec:.0f} mL")
            print("-" * 50)
            if ri >= 0.5:
                print("🟢 高可復張性 (High Recruiter)! 適合中高 PEEP 滴定，能解除塌陷並保護右心。")
            else:
                print("🚨 低可復張性 (Low Recruiter)! 強烈建議限制為中低 PEEP (8-12 cmH2O)，避免 RM 與過度膨脹。")
        except ValueError:
            print("\n[錯誤] 輸入數值格式不正確。")
            
    elif choice == '4':
        print("\n" + "="*50)
        print("          自主呼吸努力與驅力監測")
        print("="*50)
        try:
            p01 = float(input("P0.1 氣道遮斷壓 (cmH2O) [例如 2.0]: ").strip())
            dp_occ = float(input("ΔPocc 呼氣末遮斷努力壓差 (cmH2O) [例如 -12]: ").strip())
            p_peak = float(input("Peak Pressure (cmH2O) [例如 18]: ").strip())
            p_plat = float(input("Plateau Pressure (cmH2O) [例如 22]: ").strip())
            pmi = p_plat - p_peak
            
            print("-" * 50)
            print(" 安全區間對照:")
            print("  - P0.1 : 正常 < 3.5 cmH2O | 🚨 過高 >= 3.5~4.0 cmH2O")
            print("  - ΔPocc: 正常 >= -20 cmH2O | 🚨 努力過度 < -20 cmH2O")
            print("  - PMI  : 正常 <= 3.0 cmH2O | 🚨 做功過大 > 3.0 cmH2O")
            print("-" * 50)
            print(f" P0.1 狀態: {'🚨 驅力過高' if p01 >= 3.5 else '🟢 正常驅力'}")
            print(f" ΔPocc 狀態: {'🚨 努力過度 (PSILI 高危)' if dp_occ < -20.0 else '🟢 安全努力'}")
            print(f" PMI 指標: {pmi:.1f} cmH2O (PMI > 3 提示橫膈過度做功)")
            print("="*50 + "\n")
        except ValueError:
            print("\n[錯誤] 輸入數值格式不正確。")

if __name__ == "__main__":
    if "streamlit" in sys.modules or (len(sys.argv) > 1 and sys.argv[1] == "streamlit"):
        run_streamlit()
    else:
        try:
            import streamlit as st
            if st.runtime.exists():
                run_streamlit()
            else:
                run_cli()
        except ImportError:
            run_cli()
