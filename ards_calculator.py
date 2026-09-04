# -*- coding: utf-8 -*-
"""
ARDS Bedside Physiological Calculator & Dead-Space Estimator
臨床床邊生理指標與死腔預估計算器 (Streamlit Web App & CLI 雙模工具)
Based on Nuckton 2002 (NEJM) and Sinha 2019 (AJRCCM)
"""

import sys
import math

# We will implement a dual-mode script.
# If run via "streamlit run ards_calculator.py", it will execute the Streamlit app.
# If run via "python3 ards_calculator.py", it will run as a beautiful terminal interactive CLI.

# ---------------------------------------------------------
# Physiological Core Calculations
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
            "action": "1. 密切監測 PEEP 調整對 VR 的動態影響，尋找最佳 PEEP。\n2. 注意有無小氣道塌陷或早期肺部順應性變差。"
        }
    else:
        # Severe or extreme
        if vr >= 2.25:
            return {
                "level": "RED_EXTREME",
                "title": "極高度危殆區 (Extreme Risk - Quintile 5)",
                "color": "#9b0000",
                "desc": "死腔比例極高 (Vd/Vt >= 0.69，落入 Nuckton 研究最嚴重的第五個五分位區間)。",
                "action": "1. 強烈建議排除肺微血管微血栓 (Microthrombi) 或血管病變。\n2. 立即安排心臟超音波評估右心功能，防範急性肺心症 (ACP)。\n3. 考慮實施更嚴格的肺與右心保護 (Pplat < 26-28 cmH2O, pH > 7.25)。"
            }
        else:
            return {
                "level": "RED",
                "title": "高度危險區 (High Risk)",
                "color": "#e74c3c",
                "desc": "死腔比例已達或超過 ARDS 死亡組的平均值 (Vd/Vt >= 0.63)。通氣效率嚴重受損。",
                "action": "1. 評估右心負荷與肺血管阻力 (PVR)。\n2. 優化水分平衡管理 (FACCT 限制性補水/利尿)。\n3. 評估是否需啟動早期俯臥通氣 (Prone) 以優化 V/Q 匹配。"
            }

# ---------------------------------------------------------
# Streamlit App Mode
# ---------------------------------------------------------
def run_streamlit():
    import streamlit as st
    
    st.set_page_config(
        page_title="ARDS Bedside Physiological Calculator",
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
    }
    .metric-value {
        font-size: 36px;
        font-weight: bold;
        color: #2c3e50;
    }
    .metric-label {
        font-size: 14px;
        color: #7f8c8d;
        font-weight: 500;
    }
    .warning-box {
        padding: 20px;
        border-radius: 10px;
        color: white;
        margin-top: 20px;
        margin-bottom: 20px;
    }
    </style>
    """, unsafe_allow_html=True)
    
    st.markdown('<div class="main-title">🫁 ARDS 床邊生理指標與死腔預估計算器</div>', unsafe_allow_html=True)
    st.markdown('<div class="subtitle">基於 Nuckton 2002 (NEJM) 經典研究與 Sinha 2019 (AJRCCM) 通氣比例公式</div>', unsafe_allow_html=True)
    
    st.sidebar.header("📋 1. 輸入病患基本資料")
    input_method = st.sidebar.radio("PBW (預估體重) 輸入方式", ["由身高性別計算", "直接手動輸入"])
    
    pbw = 60.0
    if input_method == "由身高性別計算":
        sex = st.sidebar.selectbox("病患性別", ["男 (Male)", "女 (Female)"])
        height = st.sidebar.slider("病患身高 (cm)", 120, 210, 165)
        sex_str = "male" if "男" in sex else "female"
        pbw = calculate_pbw(height, sex_str)
        st.sidebar.info(f"計算得預估體重 (PBW): **{pbw:.1f} kg**")
    else:
        pbw = st.sidebar.number_input("預估體重 PBW (kg)", min_value=10.0, max_value=200.0, value=60.0, step=1.0)
        
    st.sidebar.header("🎛️ 2. 輸入呼吸器與血氧參數")
    ve = st.sidebar.number_input("實際分鐘通氣量 Ve (L/min)", min_value=1.0, max_value=40.0, value=12.0, step=0.5)
    paco2 = st.sidebar.number_input("實際動脈血 PaCO2 (mmHg)", min_value=10.0, max_value=150.0, value=50.0, step=1.0)
    
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
        <h3 style="margin-top:0; color:white; font-size:20px;">🚨 {risk['title']}</h3>
        <p style="margin-bottom:8px; font-size:15px; font-weight:500;">{risk['desc']}</p>
        <hr style="border-top: 1px solid rgba(255,255,255,0.3); margin: 10px 0;">
        <p style="margin-bottom:0; font-size:13px; line-height:1.5;"><strong>床邊臨床決策建議：</strong><br>{risk['action'].replace('\n', '<br>')}</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Physiological context and trial insights
    st.subheader("📚 臨床生理與研究背景 (Physiological Background)")
    
    st.markdown(r"""
    *   **死腔比例 ($V_d/V_t$) 的生死關聯**：  
        Nuckton 2002 年發表於《新英格蘭醫學期刊 (NEJM)》的經典研究證實，在 ARDS 早期，**死腔比例是強烈且完全獨立於血氧 ($PaO_2/FiO_2$) 的死亡預言家**。死腔比例每增加 **0.05**，病患的院內死亡率暴增 **45%** (Odds Ratio 1.45)。
    *   **為什麼不直接測 $V_d/V_t$ 而是算 VR？**  
        Bohr/Enghoff 經典公式需要收集 5 分鐘的混合呼出氣體量測 $P_E CO_2$，這在臨床床邊極難常規執行，且**不能直接用 $EtCO_2$ 代替 $P_E CO_2$**。而 **Ventilatory Ratio (VR)** 僅需動脈血氣的 $PaCO_2$ 與呼吸器上的分鐘通氣量 $V_E$，即可在床邊一秒算出，與真實死腔具有高度相關。
    *   **VR 數值對應標準 (Nuckton 原始數據對接)**：
        *   **VR < 1.5**：安全區，等同於 ARDS 存活者的平均死腔比例。
        *   **VR 1.63**：死亡風險開始陡增的分水嶺。
        *   **VR $\ge$ 1.90**：達到 ARDS 死亡組的平均死腔比例 (Vd/Vt >= 0.63)。
        *   **VR $\ge$ 2.26**：極高危區 (Vd/Vt >= 0.69，落入最嚴重的 20% 患者區間)。
    """)
    
    st.markdown("---")
    st.caption("聲明：本工具僅供臨床醫學學術討論與生理機制模擬使用，實際呼吸器設定與病人處置應由專科醫師依病患臨床即時狀態做出決定。")

# ---------------------------------------------------------
# CLI Command Line Mode
# ---------------------------------------------------------
def run_cli():
    print("="*65)
    print("      ARDS BEDSIDE PHYSIOLOGICAL CALCULATOR & DEAD-SPACE ESTIMATOR")
    print("          基於 Nuckton 2002 (NEJM) & Sinha 2019 生理公式")
    print("="*65)
    
    try:
        # Prompt user inputs
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
        
        # Color definition for terminal outputs (ANSI Colors)
        COLOR_MAP = {
            "GREEN": "\033[92m",
            "YELLOW": "\033[93m",
            "RED": "\033[91m",
            "RED_EXTREME": "\033[41m\033[37m" # Red background, white text
        }
        RESET = "\033[0m"
        BOLD = "\033[1m"
        
        color = COLOR_MAP.get(risk["level"], "")
        
        # Output results
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
    # If the script is called with streamlit run, let streamlit handle.
    if "streamlit" in sys.modules or (len(sys.argv) > 1 and sys.argv[1] == "streamlit"):
        run_streamlit()
    else:
        # Check if the user is calling this under a streamlit run context implicitly
        try:
            import streamlit as st
            # If we are inside streamlit run, this won't throw an error when running functions
            if st.runtime.exists():
                run_streamlit()
            else:
                run_cli()
        except ImportError:
            run_cli()
