# -*- coding: utf-8 -*-
"""
ARDS Bedside Physiological Calculator & Dead-Space Estimator v4
臨床床邊生理指標、死腔預估與生物臨床表型預測計算器 (Streamlit Web App & CLI 雙模工具)
Based on:
1. Nuckton 2002 (NEJM) & Sinha 2019 (AJRCCM) - Ventilatory Ratio (VR) & Dead Space
2. Sinha 2020 (Lancet Respiratory Medicine) - Biological Phenotype Parsimonious Models (Type 1 vs Type 2)
"""

import sys
import math

# ---------------------------------------------------------
# Physiological Core Calculations (Tab 1)
# ---------------------------------------------------------

def calculate_pbw(height_cm, sex):
    """Calculate Ideal Body Weight / Predicted Body Weight (PBW) in kg."""
    if sex.lower() in ['male', 'm', '男']:
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

# ---------------------------------------------------------
# Sinha 2020 Phenotype Parsimonious Models Calculations (Tab 2)
# ---------------------------------------------------------

# 🏆 Primary Models from main paper (Figure 2)
def calculate_primary_il8_4var_prob(il8, protein_c, bicarbonate, vasopressor_used):
    """
    Sinha 2020 Lancet RM - Primary 4-variable model (IL-8 + Protein C + Bicarbonate + Vasopressor)
    y = 4.8127 + 1.5424 * ln(IL8 + 1) - 0.2502 * Bicarbonate - 1.8778 * ln(Protein C + 1) + 2.1026 * Vasopressor
    """
    vp = 1.0 if vasopressor_used else 0.0
    log_il8 = math.log(il8 + 1.0)
    log_prot_c = math.log(protein_c + 1.0)
    y = 4.8127 + 1.5424 * log_il8 - 0.2502 * bicarbonate - 1.8778 * log_prot_c + 2.1026 * vp
    return 1.0 / (1.0 + math.exp(-y))

def calculate_primary_il8_3var_prob(il8, protein_c, bicarbonate):
    """
    Sinha 2020 Lancet RM - Primary 3-variable model (IL-8 + Protein C + Bicarbonate)
    y = 6.1241 + 1.4226 * ln(IL8 + 1) - 0.2596 * Bicarbonate - 1.8330 * ln(Protein C + 1)
    """
    log_il8 = math.log(il8 + 1.0)
    log_prot_c = math.log(protein_c + 1.0)
    y = 6.1241 + 1.4226 * log_il8 - 0.2596 * bicarbonate - 1.8330 * log_prot_c
    return 1.0 / (1.0 + math.exp(-y))

# 🧬 Ancillary Models from Supplementary Appendix Table S4
def calculate_ancillary_il6_4var_prob(il6, protein_c, bicarbonate, vasopressor_used):
    """
    Model 5 (Ancillary): IL-6 + Protein C + Bicarbonate + Vasopressor
    y = 4.0323 + 0.9191 * ln(IL6 + 1) - 0.2581 * Bicarbonate - 1.3805 * ln(Protein C + 1) + 1.7412 * Vasopressor
    """
    vp = 1.0 if vasopressor_used else 0.0
    log_il6 = math.log(il6 + 1.0)
    log_prot_c = math.log(protein_c + 1.0)
    y = 4.0323 - 0.2581 * bicarbonate + 1.7412 * vp - 1.3805 * log_prot_c + 0.9191 * log_il6
    return 1.0 / (1.0 + math.exp(-y))

def calculate_ancillary_il8_sTNFR1_4var_prob(il8, stnfr1, bicarbonate, vasopressor_used):
    """
    Model 2 (Ancillary): IL-8 + sTNFR1 + Bicarbonate + Vasopressor
    y = -13.1351 + 1.3947 * ln(IL8 + 1) - 0.2145 * Bicarbonate + 1.1818 * ln(sTNFR1 + 1) + 2.1398 * Vasopressor
    """
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
        page_title="ARDS Bedside Physiological & Phenotype Calculator",
        page_icon="🫁",
        layout="centered"
    )
    
    # Custom CSS for nicer layout and color bars
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
        font-size: 34px;
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
        padding-top: 20px;
    }
    </style>
    """, unsafe_allow_html=True)
    
    st.markdown('<div class="main-title">🫁 ARDS 床邊生理指標、死腔與生物表型計算器</div>', unsafe_allow_html=True)
    st.markdown('<div class="subtitle">整合 Nuckton 2002 (NEJM)、Sinha 2019 (AJRCCM) 與 Sinha 2020 (Lancet RM) 臨床經典研究</div>', unsafe_allow_html=True)
    
    # Define Tabs
    tab1, tab2 = st.tabs(["🫁 1. 通氣比值與死腔估算", "🧬 2. 臨床生物亞型預測"])
    
    with tab1:
        st.markdown('<div class="tab-content"></div>', unsafe_allow_html=True)
        st.subheader("📋 1. 輸入病患基本資料")
        
        # We put sidebar/main inputs
        input_method = st.radio("PBW (預估體重) 輸入方式", ["由身高性別計算", "直接手動輸入"], key="input_method")
        
        col_pbw1, col_pbw2 = st.columns(2)
        pbw = 60.0
        if input_method == "由身高性別計算":
            with col_pbw1:
                sex = st.selectbox("病患性別", ["男 (Male)", "女 (Female)"], key="sex")
            with col_pbw2:
                height = st.slider("病患身高 (cm)", 120, 210, 165, key="height")
            sex_str = "male" if "男" in sex else "female"
            pbw = calculate_pbw(height, sex_str)
            st.info(f"計算得預估體重 (PBW): **{pbw:.1f} kg**")
        else:
            pbw = st.number_input("預估體重 PBW (kg)", min_value=10.0, max_value=200.0, value=60.0, step=1.0, key="pbw_manual")
            
        st.subheader("🎛️ 2. 輸入呼吸器與血氧參數")
        col_resp1, col_resp2 = st.columns(2)
        with col_resp1:
            ve = st.number_input("實際分鐘通氣量 Ve (L/min)", min_value=1.0, max_value=40.0, value=12.0, step=0.5, key="ve")
        with col_resp2:
            paco2 = st.number_input("實際動脈血 PaCO2 (mmHg)", min_value=10.0, max_value=150.0, value=50.0, step=1.0, key="paco2")
        
        # Core calculations
        vr = calculate_vr(ve, paco2, pbw)
        vd_vt = estimate_vd_vt(vr)
        risk = get_risk_status(vr, vd_vt)
        
        # Layout of calculations
        col1, col2 = st.columns(2)
        with col1:
            st.markdown(f"""
            <div class="metric-container">
                <div class="metric-label">通氣比例 (Ventilatory Ratio, VR)</div>
                <div class="metric-value">{vr:.2f}</div>
                <div style="font-size: 11px; color:#bdc3c7;">正常健康成人為 1.0</div>
            </div>
            """, unsafe_allow_html=True)
            
        with col2:
            st.markdown(f"""
            <div class="metric-container">
                <div class="metric-label">預估生理死腔比例 (Vd/Vt)</div>
                <div class="metric-value">{vd_vt:.2f}</div>
                <div style="font-size: 11px; color:#bdc3c7;">健康上限為 0.30 | ARDS 平均 0.58</div>
            </div>
            """, unsafe_allow_html=True)
            
        # Display warning block
        bg_color = risk["color"]
        st.markdown(f"""
        <div class="warning-box" style="background-color: {bg_color};">
            <h3 style="margin-top:0; color:white; font-size:18px;">🚨 {risk['title']}</h3>
            <p style="margin-bottom:8px; font-size:14px; font-weight:500;\">{risk['desc']}</p>
            <hr style="border-top: 1px solid rgba(255,255,255,0.3); margin: 10px 0;\">
            <p style="margin-bottom:0; font-size:12px; line-height:1.5;\"><strong>床邊臨床決策建議：</strong><br>{risk['action'].replace('\\n', '<br>')}</p>
        </div>
        """, unsafe_allow_html=True)
        
        # Physiological context and trial insights
        st.subheader("📚 臨床生理與研究背景 (Physiological Background)")
        st.markdown(r"""
        *   **死腔比例 ($V_d/V_t$) 的生死關聯**：  
            Nuckton 2002 年發表於《新英格蘭醫學期刊 (NEJM)》的經典研究證使，在 ARDS 早期，**死腔比例是強烈且完全獨立於血氧 ($PaO_2/FiO_2$) 的死亡預言家**。死腔比例每增加 **0.05**，病患的院內死亡率暴增 **45%** (Odds Ratio 1.45)。
        *   **為什麼不直接測 $V_d/V_t$ 而是算 VR？**  
            Bohr/Enghoff 經典公式需要收集 5 分鐘的混合呼出氣體量測 $P_E CO_2$，這在臨床床邊極難常規執行，且**不能直接用 $EtCO_2$ 代替 $P_E CO_2$**。而 **Ventilatory Ratio (VR)** 僅需動脈血氣的 $PaCO_2$ 與呼吸器上的分鐘通氣量 $V_E$，即可在床邊一秒算出，與真實死腔具有高度相關。
        """)
        
    with tab2:
        st.markdown('<div class="tab-content"></div>', unsafe_allow_html=True)
        st.subheader("🧬 2. ARDS 臨床生物發炎亞型預測 (Sinha 2020)")
        st.markdown("""
        根據 Sinha 2020 年發表於《Lancet Respiratory Medicine》的標誌性研究，ARDS 可分為**「低發炎型 (Type 1)」**與**「高發炎型 (Type 2)」**。
        本工具已內建論文中發布的**「黃金主模型 (Primary Models)」**與**「臨床輔助模型 (Ancillary Models)」**：
        """)
        
        model_type = st.radio("選擇預測指標模型", [
            "🏆 Sinha 2020 論文主模型：IL-8 & Protein C 4變數模型 (最推薦)",
            "🏆 Sinha 2020 論文主模型：IL-8 & Protein C 3變數模型 (免收縮劑)",
            "🧬 Sinha 2020 輔助模型：IL-6 & Protein C 4變數模型 (Model 5)",
            "🧬 Sinha 2020 輔助模型：IL-8 & sTNFR-1 4變數模型 (Model 2)"
        ], key="model_selection")
        
        st.markdown("---")
        
        # Form inputs
        col_bio1, col_bio2 = st.columns(2)
        
        prob = 0.0
        
        if "論文主模型：IL-8 & Protein C 4變數" in model_type:
            with col_bio1:
                il8 = st.number_input("白介素-8 IL-8 (pg/mL)", min_value=1.0, max_value=20000.0, value=30.0, step=5.0, help="低發炎型中位數為 23，高發炎型為 133", key="il8_prim4")
                protein_c = st.number_input("蛋白質 C 活性 Protein C (% control)", min_value=1.0, max_value=250.0, value=80.0, step=5.0, help="低發炎型平均為 96.0，高發炎型為 53.5", key="pc_prim4")
            with col_bio2:
                bicarbonate_bio = st.number_input("碳酸氫根 Bicarbonate (mmol/L)", min_value=2.0, max_value=50.0, value=22.0, step=1.0, help="低發炎型平均為 23.1，高發炎型為 17.3", key="bicarb_prim4")
                vasopressor_bio = st.selectbox("是否使用血管收縮劑 (Vasopressor)", ["否 (No)", "是 (Yes)"], help="高發炎型患者 66% 需使用收縮劑，低發炎型僅 18%", key="vaso_prim4")
                
            vp_bool = True if "是" in vasopressor_bio else False
            prob = calculate_primary_il8_4var_prob(il8, protein_c, bicarbonate_bio, vp_bool)
            
        elif "論文主模型：IL-8 & Protein C 3變數" in model_type:
            with col_bio1:
                il8 = st.number_input("白介素-8 IL-8 (pg/mL)", min_value=1.0, max_value=20000.0, value=30.0, step=5.0, help="低發炎型中位數為 23，高發炎型為 133", key="il8_prim3")
                protein_c = st.number_input("蛋白質 C 活性 Protein C (% control)", min_value=1.0, max_value=250.0, value=80.0, step=5.0, help="低發炎型平均為 96.0，高發炎型為 53.5", key="pc_prim3")
            with col_bio2:
                bicarbonate_bio = st.number_input("碳酸氫根 Bicarbonate (mmol/L)", min_value=2.0, max_value=50.0, value=22.0, step=1.0, help="低發炎型平均為 23.1，高發炎型為 17.3", key="bicarb_prim3")
                st.markdown("<br><p style='font-size:12px; color:#7f8c8d; font-style:italic;'>💡 3變數模型已排除收縮劑影響，更具備跨機構決策穩定度。</p>", unsafe_allow_html=True)
                
            prob = calculate_primary_il8_3var_prob(il8, protein_c, bicarbonate_bio)
            
        elif "IL-6 & Protein C 4變數" in model_type:
            with col_bio1:
                il6 = st.number_input("白介素-6 IL-6 (pg/mL)", min_value=1.0, max_value=100000.0, value=150.0, step=10.0, help="低發炎型中位數為 116，高發炎型為 933", key="il6_anc4")
                protein_c = st.number_input("蛋白質 C 活性 Protein C (% control)", min_value=1.0, max_value=250.0, value=80.0, step=5.0, help="低發炎型平均為 96.0，高發炎型為 53.5", key="pc_anc4")
            with col_bio2:
                bicarbonate_bio = st.number_input("碳酸氫根 Bicarbonate (mmol/L)", min_value=2.0, max_value=50.0, value=22.0, step=1.0, help="低發炎型平均為 23.1，高發炎型為 17.3", key="bicarb_anc4")
                vasopressor_bio = st.selectbox("是否使用血管收縮劑 (Vasopressor)", ["否 (No)", "是 (Yes)"], key="vaso_anc4")
                
            vp_bool = True if "是" in vasopressor_bio else False
            prob = calculate_ancillary_il6_4var_prob(il6, protein_c, bicarbonate_bio, vp_bool)
            
        else: # IL-8 & sTNFR-1
            with col_bio1:
                il8 = st.number_input("白介素-8 IL-8 (pg/mL)", min_value=1.0, max_value=20000.0, value=30.0, step=5.0, help="低發炎型中位數為 23，高發炎型為 133", key="il8_anc2")
                stnfr1 = st.number_input("可溶性腫瘤壞死因子受體-1 sTNFR-1 (pg/mL)", min_value=100.0, max_value=100000.0, value=3500.0, step=100.0, help="低發炎型中位數為 3225，高發炎型為 7452", key="stnfr_anc2")
            with col_bio2:
                bicarbonate_bio = st.number_input("碳酸氫根 Bicarbonate (mmol/L)", min_value=2.0, max_value=50.0, value=22.0, step=1.0, help="低發炎型平均為 23.1，高發炎型為 17.3", key="bicarb_anc2")
                vasopressor_bio = st.selectbox("是否使用血管收縮劑 (Vasopressor)", ["否 (No)", "是 (Yes)"], key="vaso_anc2")
                
            vp_bool = True if "是" in vasopressor_bio else False
            prob = calculate_ancillary_il8_sTNFR1_4var_prob(il8, stnfr1, bicarbonate_bio, vp_bool)
            
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
                <p style="margin-bottom:8px; font-size:14px; font-weight:500;\">此患者展現出強烈的全身性發炎反應，預後風險極高 (死亡率與器官衰竭天數顯著高於 Type 1)。</p>\n                <hr style="border-top: 1px solid rgba(255,255,255,0.3); margin: 10px 0;\">
                <p style="margin-bottom:0; font-size:12px; line-height:1.5;\">
                    <strong>💡 「預測性富集」精準醫學決策價值：</strong><br>
                    1. <strong>辛伐他汀 (Simvastatin) 治療反應極佳</strong>：HARP-2 試驗事後分析證使，雖然全體分析為陰性，但 <strong>Type 2 高發炎型患者接受 Simvastatin 治療能顯著提升 28 天存活率</strong>。<br>
                    2. <strong>水分管理可能耐受寬鬆</strong>：FACCT 試驗事後分析顯示，相較於 Type 1 對保守水分有良好反應，Type 2 患者在特定情況下對寬鬆補液可能呈現更好的生存趨勢。<br>
                    3. <strong>高度右心衰竭風險</strong>：伴隨嚴重酸縮與代謝不全，應極度嚴防急性肺心症 (ACP)。
                </p>
            </div>
            """, unsafe_allow_html=True)
        else:
            st.markdown(f"""
            <div class="warning-box" style="background-color: #2ebd59;">
                <h3 style="margin-top:0; color:white; font-size:18px;">🧬 預測亞型：Type 1 低發炎型 (Hypoinflammatory)</h3>
                <p style="margin-bottom:8px; font-size:14px; font-weight:500;\">全身性發炎反應較為輕微，預後顯著優於高發炎型 (約佔全體 ARDS 患者的 70%)。</p>\n                <hr style="border-top: 1px solid rgba(255,255,255,0.3); margin: 10px 0;\">
                <p style="margin-bottom:0; font-size:12px; line-height:1.5;\">
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
            研究最初使用包含 30 餘種臨床和分子指標的 LCA 分析。為了實用化，Sinha 2020 年開發了這套「簡化模型」——透過 3-4 個關鍵指標（IL-8, Bicarbonate, Protein C, Vasopressor），可在不損失預測力 (AUC 達 0.94 - 0.96) 的前提下在床邊快速分類，對未來 ARDS 的精準用藥 (Targeted Therapy) 具有奠基地位。
        """)

    st.markdown("---")
    st.caption("聲明：本工具僅供臨床醫學學術討論與生理機制模擬使用，實際呼吸器設定與病人處置應由專科醫師依病患臨床即時狀態做出決定。")

# ---------------------------------------------------------
# CLI Command Line Mode
# ---------------------------------------------------------
def run_cli():
    print("="*65)
    print("      ARDS BEDSIDE PHYSIOLOGICAL & BIOLOGICAL PHENOTYPE CALCULATOR")
    print("                臨床生理、死腔與生物臨床表型計算器 v4")
    print("="*65)
    
    print("\n請選擇您要執行的功能：")
    print(" [1] 計算 Ventilatory Ratio (VR) 與生理死腔 (Nuckton 2002 / Sinha 2019)")
    print(" [2] 預測 ARDS 生物臨床發炎亞型 (Sinha 2020)")
    
    choice = input("\n請輸入 1 或 2: ").strip()
    
    if choice == '2':
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
            
            if m_choice in ['1', '3', '4']:
                vaso_input = input("是否使用血管收縮劑? (y/n): ").strip().lower()
                vp_bool = True if vaso_input in ['y', 'yes', '是'] else False
            else:
                vp_bool = False
            
            if m_choice == '1':
                il8 = float(input("請輸入 IL-8 (pg/mL) [例如 30]: ").strip())
                protein_c = float(input("請輸入 Protein C (% control) [例如 80]: ").strip())
                prob = calculate_primary_il8_4var_prob(il8, protein_c, bicarbonate, vp_bool)
                model_name = "Sinha 2020 主模型 (IL-8 & Protein C 4變數)"
            elif m_choice == '2':
                il8 = float(input("請輸入 IL-8 (pg/mL) [例如 30]: ").strip())
                protein_c = float(input("請輸入 Protein C (% control) [例如 80]: ").strip())
                prob = calculate_primary_il8_3var_prob(il8, protein_c, bicarbonate)
                model_name = "Sinha 2020 主模型 (IL-8 & Protein C 3變數)"
            elif m_choice == '3':
                il6 = float(input("請輸入 IL-6 (pg/mL) [例如 150]: ").strip())
                protein_c = float(input("請輸入 Protein C (% control) [例如 80]: ").strip())
                prob = calculate_ancillary_il6_4var_prob(il6, protein_c, bicarbonate, vp_bool)
                model_name = "Sinha 2020 輔助模型 (IL-6 & Protein C 4變數)"
            else:
                il8 = float(input("請輸入 IL-8 (pg/mL) [例如 30]: ").strip())
                stnfr1 = float(input("請輸入 sTNFR-1 (pg/mL) [例如 3500]: ").strip())
                prob = calculate_ancillary_il8_sTNFR1_4var_prob(il8, stnfr1, bicarbonate, vp_bool)
                model_name = "Sinha 2020 輔助模型 (IL-8 & sTNFR-1 4變數)"
                
            prob_pct = prob * 100.0
            
            # ANSI Colors
            RED = "\033[91m"
            GREEN = "\033[92m"
            RESET = "\033[0m"
            BOLD = "\033[1m"
            
            print("\n" + "="*50)
            print(f"{BOLD}預測結果 (Prediction Results):{RESET}")
            print(f" 使用模型: {model_name}")
            print(f" 高發炎亞型 (Type 2) 機率: {BOLD}{prob_pct:.1f}%{RESET}")
            print("-"*50)
            
            if prob >= 0.5:
                print(f"🚨 {RED}{BOLD}預測亞型：Type 2 高發炎型 (Hyper-inflammatory){RESET}")
                print(" 生理狀態: 全身發炎反應強烈，預後風險極高。")
                print(" 臨床決策價值:\n  1. HARP-2 分析證實 Simvastatin 能顯著提升 Type 2 存活率。\n  2. 對水分寬鬆補液策略可能呈現較佳反應。")
            else:
                print(f"🟢 {GREEN}{BOLD}預測亞型：Type 1 低發炎型 (Hypoinflammatory){RESET}")
                print(" 生理狀態: 發炎反應較輕，預後顯著優異。")
                print(" 臨床決策價值:\n  1. 保守水分限制策略 (FACCT 策略) 能最大幅減少呼吸器天數。\n  2. Simvastatin 治療對此類病患無益，甚至可能有副作用。")
            print("="*50 + "\n")
            
        except ValueError:
            print("\n[錯誤] 輸入數值格式不正確，請重新執行。")
            
    else:
        # Standard VR calculation CLI
        try:
            print("\n[1] 病患預估體重 (PBW) 計算")
            calc_choice = input("是否需要從身高與性別計算 PBW? (y/n, 預設為直接手動輸入): ").strip().lower()
            
            if calc_choice == 'y':
                sex_input = input("病患性別 (m=男, f=女): ").strip().lower()
                sex = 'male' if sex_input in ['m', 'male', '男'] else 'female'
                height = float(input("病患身高 (cm): ").strip())
                pbw = calculate_pbw(height, sex)
                print(f"--> 計算得出預估體重 (PBW): {pbw:.2f} kg")
            else:
                pbw = float(input("請輸入預估體重 PBW (kg): ").strip())
                
            print("\n[2] 呼吸器與血氧參數")
            ve = float(input("實際分鐘通氣量 Ve (L/min): ").strip())
            paco2 = float(input("實際動脈血 PaCO2 (mmHg): ").strip())
            
            # Calculations
            vr = calculate_vr(ve, paco2, pbw)
            vd_vt = estimate_vd_vt(vr)
            risk = get_risk_status(vr, vd_vt)
            
            COLOR_MAP = {
                "GREEN": "\033[92m",
                "YELLOW": "\033[93m",
                "RED": "\033[91m",
                "RED_EXTREME": "\033[41m\033[37m"
            }
            RESET = "\033[0m"
            BOLD = "\033[1m"
            
            color = COLOR_MAP.get(risk["level"], "")
            
            print("\n" + "="*50)
            print(f"{BOLD}計算結果 (Bedside Results):{RESET}")
            print(f" - 通氣比例 (Ventilatory Ratio, VR): {BOLD}{vr:.2f}{RESET} (健康成人正常值為 1.0)")
            print(f" - 預估生理死腔比例 (Vd/Vt Estimator): {BOLD}{vd_vt:.2f}{RESET} (健康上限為 0.30)")
            print("-"*50)
            print(f"🚨 警訊分級: {color}{BOLD}{risk['title']}{RESET}")
            print(f"📊 生理狀態: {risk['desc']}")
            print(f"👉 臨床決策建議:\n{risk['action']}")
            print("="*50 + "\n")
            
        except ValueError:
            print("\n[錯誤] 輸入數值格式不正確，請重新執行並輸入數字。")
        except KeyboardInterrupt:
            print("\n\n計算器已關閉。")

# ---------------------------------------------------------
# Entry Point
# ---------------------------------------------------------
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
