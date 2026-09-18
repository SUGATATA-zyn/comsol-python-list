# -*- coding: utf-8 -*-
"""
WSL 闭环外翻滚动 · 步骤 2 规范算例:轴对称充压段打靶  (shoot_free_arc.py, v7)
================================================================================

几何(v7 修正,依据《轴对称封闭膜结构——几何理解与建模指南》):
  * 该结构是 **二维 (r,z) 平面中、绕 r=0 轴对称的"封闭膜带"**: 膜不是单条中
    心线, 而是由 **膜外侧轮廓** 与 **膜内侧轮廓** 两条闭合曲线夹成的 1 mm 厚
    膜实体; 膜内侧轮廓所围区域 = **完全封闭的内部腔体**(不是外部流体域)。
  * 左竖直膜段(靠近对称轴一侧): 膜外侧 r=14 mm、膜内侧 r=15 mm(即外侧比
    内侧更靠近 r=0; 常见 "R_out>R_in" 的理解在这里是错的)。左侧直线段长
    H1=90 mm, 右侧直线段长 H2=120 mm, 膜厚 t=1 mm。
  * 上端过渡(从左壁到右壁, 逐段相切、至少 G1): 左侧竖直直线 -> 半径 15 mm
    的 45° 圆弧 -> 45° 斜直线(L=21.213 mm, dr=dz=15 mm) -> 半径 15 mm 的
    135° 圆弧 -> 右侧竖直直线; 下端为同一三段式的反向(镜像)连接。
  * 右侧外轮廓与独立外部矩形计算域(刚性约束壁)左边界相距 1 mm; 该刚性壁即
    约束半径 R_w=62 mm。充压后膜外侧右壁(r≈61 mm)压贴到 r=62 mm 的壁面,
    贴壁段由此形成(摩擦/CZM 切向锚定), 滚轮驱动使前端从壁面剥离并内翻。
  * 降阶膜理论取 **中面(距内外轮廓各 0.5 mm)闭合成环** 作参考构型, 以参考
    弧长 S in [0, L0] 参数化 r0(S)、z0(S)。中面周长由几何参数唯一确定:
        L0 = H1 + H2 + 2*21.213 + 2*pi*R_mid   (= 349.815 mm)
    其中 R_mid=15.5 mm(内弧 R=15、外弧 R=16 的中间面, 圆弧与内弧同圆心), 与
    几何模型/wsl_membrane_geom.py 生成并自检通过的轮廓一致(左 14/15、右 60/61)。
  * 参考环分段表(单位 mm, 坐标为物理 r/z, 关于 z=0 对称; S=0 暂取左直壁下端点
    (r,z)=(14.5,-45) 逆时针成环; 与 ODE“贴壁段从 S=0 开始”的拼接见下“待标定”):
        [0, 90]            左侧竖直线(上行)     r=14.5,   z: -45 -> 45
        [90, ~102.17]      上端 45° 圆弧 R15.5  r: 14.5 -> 19.04, z: 45 -> 55.96
        [~102.17, ~123.39] 45° 斜线 21.213     dr=dz=+15 (r: 19.04 -> 34.04)
        [~123.39, ~159.91] 上端 135° 圆弧 R15.5 r: 34.04 -> 60.5 (过顶 z=75.5)
        [~159.91, ~279.91] 右侧竖直线(下行)    r=60.5,   z: 60 -> -60
        [~279.91, ~316.43] 下端 135° 圆弧 R15.5 r: 60.5 -> 34.04 (过底 z=-75.5)
        [~316.43, ~337.64] 45° 斜线 21.213     dr=-15, dz=+15 (r: 34.04 -> 19.04)
        [~337.64, L0]      下端 45° 圆弧 R15.5 r: 19.04 -> 14.5 (回起点 z=-45)
    圆弧圆心与内侧轮廓同圆心: 45° 弧 (r,z)=(30,±45)、135° 弧 (45,±60) mm,
    中面弧半径 R_mid=15.5 mm; 45° 投影精确 dr=dz=15 mm; 弧长=R_mid*转角。
    r0(S)/z0(S) 由 ref_r0/ref_z0 逐段解析给出(圆弧解析, 非按弧长线性插值)。

本构:不可压缩 Yeoh(平面应力膜),张力-主伸长(推导笔记式 2.5):
      N_a = 2 h0 W1(I1) lam3 (lam_a^2 - lam3^2),  lam3 = 1/(lam1 lam2)
自由段平衡(以参考弧长 S 积分,内压 p 沿腔侧外法向;符号约定已在 §3.1 校验):
      dr/dS   = lam1 cos(psi)
      dz/dS   = lam1 sin(psi)
      dpsi/dS = lam1*(+p - N_phi*sin(psi)/r)/N_s
      dN_s/dS = -lam1*(N_s - N_phi)*cos(psi)/r
贴壁段: 膜中面压贴壁面(刚性壁半径 R_w=62 mm -> 中面接触半径 R_c=R_w-t/2,
  切线 psi=pi/2); N_s 常数; 接触压力 q(S)=p+N_phi/R_c。

关于定解(重要,对应推导笔记 §3.4 的展开):
  纯膜 + 无黏附 + 切线连续剥离是**过定**的: 两个未知量 {N_s0, s_c} 要同时满足
  闭端 z=0、psi=-pi/2 与剥离点 q=0 三个条件; 缺的自由度正是"剥离角/剥离几何"。
  COMSOL 中贴壁段切向锚定(摩擦+CZM)与剥离处弯曲边界层补足该自由度, 使 s_c 成
  为由加载历史决定的量。本代码把 **s_c(贴壁段材料长度)作为输入/扫描参数**,
  未知量 x=(N_s0, psi_c), 条件 z(S=L0/2)=0 与
  wrap_to_pi(psi(S=L0/2)-3*pi/2)=0。
  COMSOL para=1 材料坐标初步标定: 全局参考环 S=0 保持在左壁下端, 但半环 ODE
  使用独立的局部材料坐标 S_hat。S_hat=0 是右壁中面 z=0 的材料点, 并沿上半环
  方向增加; 因此贴壁段 [0,s_c] 与其真实参考材料区间一致。局部到全局映射由
  _SEG 中 tag='R' 的右壁段自动导出, 不改变全局参考几何。

  数值诊断: 在 (Ns0, psi_c) 平面上残差 r=(z_end, angle_res) 呈**近平行双峡谷**
  结构——z_end=0 与 angle_res=0 两条零曲线在相关参数域内接近但通常不相交, 因此
  对给定 (p, s_c) **精确纯膜闭解往往不存在**, 只存在"最小残差准闭解"(峡谷
  底)。峡谷间距即"纯膜定解缺掉的那个自由度"的数值体现: 真实结构中由弯曲边界
  层/摩擦锚定/CZM 闭合。v7 求解器因此输出两类结果并**如实标注残差**:
    * 精确闭解: |z_end|,|angle_res| 均 ~1e-9 (rare, 需 (p,s_c) 恰在交点);
    * 准闭解(默认): 峡谷底最优点, 逐点打印残差, 供与 COMSOL 构型对照。

单位:SI。运行示例:
  python shoot_free_arc.py                          # 默认演示: p=400 Pa, s_c 扫描
  python shoot_free_arc.py --p 400 --scan-sc 8,10.5,13,16,20,25e-3
  python shoot_free_arc.py --p 400 --fix-sc 0.0105  # 单点(对照 COMSOL para=1)
  python shoot_free_arc.py --scan 300,350,400,450,500 --fix-sc 0.0105
输出:控制台表;--out 指定 CSV(逐算例写构型)。
"""
import math, argparse, csv

# ---------------- 输入参数(依据《轴对称封闭膜结构——几何理解与建模指南》) ----------------
C10 = 110000.0     # Pa   Yeoh 一阶
C20 = 20000.0      # Pa   Yeoh 二阶
C30 = 0.0          # Pa   Yeoh 三阶
H0  = 1.0e-3       # m    初始膜厚 t=1 mm(膜内外轮廓间距)

# ---- 结构尺寸(几何指南 §3-§8) ----
RM_OUT_L = 14.0e-3   # m  左侧膜段: 膜外侧半径(靠 r=0 一侧)
RM_IN_L  = 15.0e-3   # m  左侧膜段: 膜内侧半径(靠腔体一侧)
R_MID_L  = 0.5*(RM_OUT_L + RM_IN_L)   # 左竖直壁中面半径 = 14.5 mm
H1       = 90.0e-3   # m  左侧竖直直线段长度
H2       = 120.0e-3  # m  右侧竖直直线段长度
ARC_R    = 15.0e-3   # m  圆弧半径 R=15 mm(指主/内侧轮廓圆弧)
R_MID_A  = ARC_R + H0/2.0             # 中面圆弧半径 = 15.5 mm(闭合周长校验用)
D_45     = 21.213e-3 # m  45° 斜线长度(L=21.213 mm)
D_RZ     = D_45/math.sqrt(2.0)       # 斜线在 r/z 向投影 = 15 mm

# ---- 刚性约束壁(几何指南 §11: 右侧矩形域左边界) ----
RW = 62.0e-3        # m  刚性圆柱壁(约束壁)半径, 膜外右侧轮廓初始距壁 1 mm
RC = RW - H0/2.0    # m  贴壁时膜中面接触半径 = 61.5 mm

# ---- COMSOL para=1 对照标定值(仅记录当前标定, 不改变 CLI 通用扫描默认值) ----
P_COMSOL_REF = 300.0       # Pa
# 原始 COMSOL 分离点:
#   upper peel: 参考 (R,Z)=(61.0,+50.85) mm, 变形后 (r,z)=(62.0,+50.685) mm;
#   lower peel: 参考 (R,Z)=(61.0,-51.00) mm, 变形后 (r,z)=(61.999,-50.835) mm.
# raw material contact length = 101.85 mm; 两端 0.15 mm 差异暂视作网格/端点识别
# 的小非对称。轴对称半环理论取 z=0 对称化后 s_c=50.925 mm。
LC_MATERIAL_REF = 0.10185  # m, 参考构型材料接触总长
LC_CURRENT_COMSOL = 0.1017 # m, 当前变形构型总接触长(不用于 s_c 定义)
SC_COMSOL_REF = 0.050925   # m, 对称化的半环参考材料贴壁长度

# ---- 中面参考环(闭合分段几何, 与 几何模型/wsl_membrane_geom.py 严格一致) ----
# 主/内轮廓 = 腔侧, 圆弧 R=ARC_R=15 mm(圆心见下); 外轮廓 = 内轮廓沿离腔侧法线偏
# 1 mm(同圆心 R=16); 中面 = 内外两侧各偏 0.5 mm 的闭合成环: 圆弧半径 R_MID_A=15.5
# 且圆心与内弧相同。圆心(按几何脚本): 45° 弧 (CR_A45,±ZC_L45)=(30,±45) mm,
# 135° 弧 (CR_A135,±ZC_R135)=(45,±60) mm。直壁: 左 r=14.5、z∈[-45,45]
# (长 90 mm), 右 r=60.5、z∈[-60,60](长 120 mm)。
CR_A45  = RM_IN_L + ARC_R      # m   45° 圆弧圆心半径 = 30 mm
CR_A135 = CR_A45 + D_RZ        # m   135° 圆弧圆心半径 = 45 mm
ZC_L45  = H1/2.0               # m   45° 弧圆心 |z| = 45 mm(左直壁半高)
ZC_R135 = H2/2.0               # m   135° 弧圆心 |z| = 60 mm(右直壁半高)
COS45   = math.cos(math.pi/4.0)
DR_A45  = R_MID_A*(1.0 - COS45)   # m  45° 弧径向增量 = 4.5398 mm
DR_A135 = R_MID_A*(1.0 + COS45)   # m  135° 弧径向增量 = 26.4602 mm
R_UP45  = R_MID_L + DR_A45        # m  45° 弧终点(=斜线起点) r = 19.0398 mm
R_UP    = R_UP45 + D_RZ           # m  斜线终点(=135° 弧起点) r = 34.0398 mm
R_MID_R = R_UP + DR_A135          # m  右直壁中面半径 = 60.5 mm(与几何脚本一致)
APEX_Z  = ZC_R135 + R_MID_A       # m  上端最高点 z = 75.5 mm(内 75 / 外 76 的中面)
Z_A45   = ZC_L45 + R_MID_A*COS45  # m  45° 弧终点 z = 55.9602 mm
Z_UP    = Z_A45 + D_RZ            # m  斜线终点 z = 70.9602 mm

def _build_segments():
    """按参考弧长 S 生成中面闭合参考环分段表(单位 m)。S=0 暂取左直壁下端
    (r,z)=(14.5,-45) mm, 逆时针成环: 左壁上行 -> 上端 45° 弧 -> 斜线 -> 上端
    135° 弧 -> 右壁下行 -> 下端 135° 弧 -> 斜线 -> 45° 弧 -> 回到起点。
    每段 dict(kind, tag, s0, s1, 几何量); 圆弧参数: 圆心 (cr,cz)、半径 R、起点角
    th0, 且 θ 随 S 递减: θ(S) = th0 - (S-s0)/R。弧长 = R_MID_A*转角,
    周长 L0 自动为 349.815 mm。"""
    a45 = math.pi/4.0
    s = 0.0
    seg = []
    seg.append(dict(kind='wall',  tag='L', s0=s, s1=s+H1, r=R_MID_L))
    s += H1
    seg.append(dict(kind='arc', tag='a45', s0=s, s1=s+R_MID_A*a45,
                    cr=CR_A45, cz=ZC_L45, R=R_MID_A, th0=math.pi))
    s += R_MID_A*a45
    seg.append(dict(kind='slant', tag='U', s0=s, s1=s+D_45,
                    r0=R_UP45, r1=R_UP, z0=Z_A45, z1=Z_UP))
    s += D_45
    seg.append(dict(kind='arc', tag='A135', s0=s, s1=s+R_MID_A*3.0*a45,
                    cr=CR_A135, cz=ZC_R135, R=R_MID_A, th0=3.0*a45))
    s += R_MID_A*3.0*a45
    seg.append(dict(kind='wall',  tag='R', s0=s, s1=s+H2, r=R_MID_R))
    s += H2
    seg.append(dict(kind='arc', tag='a135', s0=s, s1=s+R_MID_A*3.0*a45,
                    cr=CR_A135, cz=-ZC_R135, R=R_MID_A, th0=0.0))
    s += R_MID_A*3.0*a45
    seg.append(dict(kind='slant', tag='D', s0=s, s1=s+D_45,
                    r0=R_UP, r1=R_UP45, z0=-Z_UP, z1=-Z_A45))
    s += D_45
    seg.append(dict(kind='arc', tag='A45', s0=s, s1=s+R_MID_A*a45,
                    cr=CR_A45, cz=-ZC_L45, R=R_MID_A, th0=-3.0*a45))
    s += R_MID_A*a45
    return seg

_SEG = _build_segments()
L0   = _SEG[-1]['s1']            # m  中面参考总弧长 = 349.815 mm(自检)

# 半环局部材料坐标 S_hat: 0 位于右壁中面 z=0, 沿上半环方向增加。
# 全局 S 在右壁上沿 z=+60 -> -60 mm 递增, 所以局部正方向与全局方向相反。
_R_WALL_SEG = next(g for g in _SEG if g['kind'] == 'wall' and g['tag'] == 'R')
S_R_START = _R_WALL_SEG['s0']
S_CONTACT_CENTER_GLOBAL = S_R_START + 0.5*(_R_WALL_SEG['s1'] - _R_WALL_SEG['s0'])

def local_to_global_S(S_hat):
    """将局部半环材料坐标 S_hat 映射到全局参考环 S(周期处理)。"""
    return (S_CONTACT_CENTER_GLOBAL - S_hat) % L0

def ref_r0(S):
    """参考(未变形)构型中材料点 S 的中面半径 r0(S)。S∈[0,L0]; 圆弧用解析式
    r0 = cr + R*cos(th0-(S-s0)/R)(非按弧长线性插值), 保证环上 r0(S) 连续。"""
    if S < 0.0 or S > L0:
        raise ArithmeticError('ref_r0: S=%g 超出 [0, L0]' % S)
    for g in _SEG:
        if S <= g['s1']:
            if g['kind'] == 'wall':
                return g['r']
            if g['kind'] == 'slant':
                u = (S - g['s0'])/(g['s1'] - g['s0'])
                return g['r0'] + (g['r1'] - g['r0'])*u
            th = g['th0'] - (S - g['s0'])/g['R']
            return g['cr'] + g['R']*math.cos(th)
    raise ArithmeticError('ref_r0: S=%g 未命中分段' % S)

def ref_z0(S):
    """参考构型中材料点 S 的轴向坐标 z0(S)(物理坐标, 关于 z=0 对称)。
    供构型画图/对照与几何自检; ODE 只用到 r0(S)(经 lambda2 = r/r0)。"""
    if S < 0.0 or S > L0:
        raise ArithmeticError('ref_z0: S=%g 超出 [0, L0]' % S)
    for g in _SEG:
        if S <= g['s1']:
            if g['kind'] == 'wall':
                if g['tag'] == 'L':
                    return -ZC_L45 + (S - g['s0'])
                return ZC_R135 - (S - g['s0'])
            if g['kind'] == 'slant':
                u = (S - g['s0'])/(g['s1'] - g['s0'])
                return g['z0'] + (g['z1'] - g['z0'])*u
            th = g['th0'] - (S - g['s0'])/g['R']
            return g['cz'] + g['R']*math.sin(th)
    raise ArithmeticError('ref_z0: S=%g 未命中分段' % S)

def ref_r0_local(S_hat):
    """半环 ODE 的局部材料点参考中面半径 r0(S_hat)。"""
    return ref_r0(local_to_global_S(S_hat))

def ref_z0_local(S_hat):
    """半环 ODE 的局部材料点参考轴向坐标 z0(S_hat)。"""
    return ref_z0(local_to_global_S(S_hat))

def local_material_mapping_check(tol=1.0e-10):
    """校验 COMSOL para=1 对称化局部材料坐标与全局参考环的对应关系。"""
    points = {
        'center': (0.0, R_MID_R, 0.0),
        'peel': (SC_COMSOL_REF, R_MID_R, SC_COMSOL_REF),
        'half_end': (L0/2.0, R_MID_L, 0.0),
    }
    values = {}
    for name, (S_hat, expected_r, expected_z) in points.items():
        r0 = ref_r0_local(S_hat)
        z0 = ref_z0_local(S_hat)
        if abs(r0 - expected_r) > tol or abs(z0 - expected_z) > tol:
            raise AssertionError(
                'local material mapping %s failed: r0=%g z0=%g' % (name, r0, z0))
        values[name] = (r0, z0)
    return values

def geometry_report():
    """几何自检: 打印中面参考环分段端点与周长(理论 L0 应 = 349.815 mm)。
    过渡段半径应连续: 14.5 -> 19.04 -> 34.04 -> 60.5(镜像回 14.5)。"""
    Lh = _SEG[-1]['s1']
    apex = _SEG[3]['s0'] + _SEG[3]['R']*math.pi/4.0  # 上端 135° 弧最高点(S)
    lines = ['== 参考几何自检(中面环) ==',
             '  直壁中面 r: 左 = %.3f mm | 右 = %.3f mm(z: 左 ±45 / 右 ±60)'
             % (R_MID_L*1e3, R_MID_R*1e3),
             '  过渡段 r: 45° 弧后 %.4f mm -> 斜线后 %.4f mm -> 135° 弧后 %.4f mm'
             % (R_UP45*1e3, R_UP*1e3, R_MID_R*1e3),
             '  中面弧 R = %.3f mm (45+135+135+45 deg), 45 斜线 %.4f mm x2, L0 = %.4f mm'
             % (R_MID_A*1e3, D_45*1e3, Lh*1e3),
             '  最高点 z = %.4f mm (r = %.4f mm); 最低点 z = %.4f mm'
             % (ref_z0(apex)*1e3, ref_r0(apex)*1e3, -ref_z0(apex)*1e3)]
    return '\n'.join(lines)

# ============================================================
#  1. 本构:不可压缩 Yeoh,平面应力(式 2.5)
# ============================================================
def W1_of(I1):
    m = I1 - 3.0
    return C10 + 2.0*C20*m + 3.0*C30*m*m

def yeoh_N(lam1, lam2):
    """(N1,N2) 子午向/环向膜合力 [N/m]。"""
    lam3 = 1.0/(lam1*lam2)
    I1 = lam1*lam1 + lam2*lam2 + lam3*lam3
    w1 = W1_of(I1)
    c = 2.0*H0*w1*lam3
    return c*(lam1*lam1 - lam3*lam3), c*(lam2*lam2 - lam3*lam3)

def inv_lam1(N1t, lam2):
    """给 (目标张力 N1t>=0, 环向伸长 lam2>0) 反解 lam1(牛顿+回退二分)。"""
    if N1t < 0.0 or not (0.05 <= lam2 <= 8.0):
        raise ArithmeticError('inv_lam1: N1t=%g lam2=%g 越界' % (N1t, lam2))
    z = lam2**(-1.0/3.0)
    if N1t <= 0.0:
        return z
    lam = z*1.35
    for _ in range(14):
        N1v, _ = yeoh_N(lam, lam2)
        if abs(N1v - N1t) < 1e-9*max(N1t, 1.0) + 1e-12:
            return lam
        d = 1e-4*max(lam, 1e-3)
        N1d, _ = yeoh_N(lam + d, lam2)
        slope = (N1d - N1v)/d
        if slope <= 0.0 or not math.isfinite(slope):
            break
        lam -= (N1v - N1t)/slope
        if lam < 0.7*z:
            lam = 0.5*(lam + z)
        if lam > 1e5:
            raise ArithmeticError('inv_lam1: N1t=%g 超出可反解范围' % N1t)
    lo = max(0.5*z, 0.02); hi = z
    for _ in range(200):
        N1h, _ = yeoh_N(hi, lam2)
        if N1h >= N1t:
            break
        hi *= 1.6
        if hi > 1e5:
            raise ArithmeticError('inv_lam1: N1t=%g 超出可反解范围' % N1t)
    for _ in range(60):
        mid = 0.5*(lo+hi)
        N1m, _ = yeoh_N(mid, lam2)
        if N1m < N1t:
            lo = mid
        else:
            hi = mid
    return 0.5*(lo+hi)

# ============================================================
#  2. 半环积分:贴壁段[0,s_c] + 自由段[s_c,L0/2]
# ============================================================
def wrap_to_pi(angle):
    """将角度映射到 [-pi, pi)，用于有方向的终端切线残差。"""
    return (angle + math.pi) % (2.0*math.pi) - math.pi

def integrate_half(Ns0, s_c, psi_c, p, nstep):
    """贴壁段 psi=pi/2、N_s=Ns0、r=RC;自由段起点 psi=psi_c、N_s=Ns0/sin(psi_c)。
    返回轨迹字典;自由段张力非正或越界时抛 ArithmeticError。"""
    Lh = L0/2.0
    sp = math.sin(psi_c)
    if Ns0 <= 1e-9 or not (0.0 <= s_c < Lh*0.999) or sp <= 1e-6:
        raise ArithmeticError('integrate_half: 初值不合法')
    # psi_c 当前搜索范围在 (pi/2, pi)，故 sin(psi_c)>0；该张力转换的物理
    # 推导在后续阶段单独复核，本轮保持原式。
    Nsp = Ns0/sp
    h = Lh/nstep
    sc_steps = int(round(s_c/h))
    l20 = RC/ref_r0_local(0.0)
    l10 = inv_lam1(Ns0, l20)
    _, Np20 = yeoh_N(l10, l20)
    S=[0.0]; r=[RC]; z=[0.0]; psi=[math.pi/2.0]; Ns=[Ns0]
    lam1=[l10]; lam2=[l20]; q=[p + Np20/RC]
    for i in range(sc_steps):
        Sp = (i+1)*h
        l2 = RC/ref_r0_local(Sp)
        l1 = inv_lam1(Ns0, l2)
        _, Np2 = yeoh_N(l1, l2)
        S.append(Sp); r.append(RC); z.append(z[-1]+l1*h)
        psi.append(math.pi/2.0); Ns.append(Ns0)
        lam1.append(l1); lam2.append(l2); q.append(p + Np2/RC)

    def deriv(Sv, y):
        rr, zz, ps, N1 = y
        if N1 <= 1e-9 or rr <= 0.0:
            raise ArithmeticError('deriv: Ns<=0 或 r<=0')
        l2 = rr/ref_r0_local(Sv)
        l1 = inv_lam1(N1, l2)
        _, N2 = yeoh_N(l1, l2)
        cp = math.cos(ps); sn = math.sin(ps)
        return (l1*cp, l1*sn, l1*(p - N2*sn/rr)/N1, -l1*(N1 - N2)*cp/rr)

    def rk4(Sv, y):
        k1 = deriv(Sv, y)
        k2 = deriv(Sv+h/2, [y[i]+0.5*h*k1[i] for i in range(4)])
        k3 = deriv(Sv+h/2, [y[i]+0.5*h*k2[i] for i in range(4)])
        k4 = deriv(Sv+h,   [y[i]+h*k3[i] for i in range(4)])
        return [y[i]+h*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0 for i in range(4)]

    y = [RC, z[-1], psi_c, Nsp]
    for i in range(sc_steps, nstep):
        y = rk4(i*h, y)
        Sv = (i+1)*h
        rr, zz, ps, N1 = y
        l2 = rr/ref_r0_local(Sv)
        l1 = inv_lam1(N1, l2)
        S.append(Sv); r.append(rr); z.append(zz); psi.append(ps)
        Ns.append(N1); lam1.append(l1); lam2.append(l2); q.append(0.0)
    psi_end = psi[-1]
    angle_res = wrap_to_pi(psi_end - 3.0*math.pi/2.0)
    l2_free_start = RC/ref_r0_local(sc_steps*h)
    l1_free_start = inv_lam1(Nsp, l2_free_start)
    return dict(S=S, r=r, z=z, psi=psi, Ns=Ns, lam1=lam1, lam2=lam2,
                q=q, sc_steps=sc_steps, sc_eff=sc_steps*h,
                res=(z[-1], angle_res), angle_res=angle_res,
                psi_start=psi[0], psi_c=psi_c, psi_end=psi_end,
                r_start=r[0], r_peel=r[sc_steps], r_end=r[-1],
                drds_free_start=l1_free_start*math.cos(psi_c))

# ============================================================
#  3. 求解器:全域粗采样 -> 多分支选优 -> 逐分支放大抛光
#     (未知 Ns0, psi_c; 残差 r = (z_end, angle_res))
# ============================================================
def _safe_res(Ns0, psi_c, s_c, p, nstep, cache):
    key = (round(Ns0, 11), round(psi_c, 11))
    v = cache.get(key)
    if v is None:
        try:
            v = integrate_half(Ns0, s_c, psi_c, p, nstep)['res']
        except Exception:
            v = None
        cache[key] = v
    return v

def _nrm(r):
    return abs(r[0]) + abs(r[1]) if r is not None else float('inf')

def _newton(rr, x0, maxit=80):
    """阻尼牛顿;收敛到 |r|<1e-12 返回 (x,r),否则 None。"""
    x = [float(x0[0]), float(x0[1])]
    r0 = rr(x[0], x[1])
    if r0 is None:
        return None
    for _ in range(maxit):
        if _nrm(r0) < 1e-13:
            return (x, r0)
        J = []; bad = False
        for k in range(2):
            hx = 1e-6*max(abs(x[0]), 1.0) if k == 0 else 2e-7
            xp = list(x); xp[k] += hx
            rp = rr(xp[0], xp[1])
            if rp is None:
                bad = True; break
            J.append([(rp[0]-r0[0])/hx, (rp[1]-r0[1])/hx])
        if bad:
            return None
        det = J[0][0]*J[1][1] - J[0][1]*J[1][0]
        if abs(det) < 1e-28:
            return None
        dx = [-(r0[0]*J[1][1]-r0[1]*J[0][1])/det,
              -(J[0][0]*r0[1]-J[1][0]*r0[0])/det]
        b0 = _nrm(r0)
        t = 1.0
        while t > 1e-12:
            xn = [x[0]+t*dx[0], x[1]+t*dx[1]]
            rn = rr(xn[0], xn[1])
            if rn is not None and _nrm(rn) < b0:
                x, r0 = xn, rn
                break
            t *= 0.5
        else:
            return None
    return (x, r0) if _nrm(r0) < 1e-12 else None

def _pattern_min(rr, x0, n_step0=(0.02, 0.012), maxit=240):
    """交替方向模式下降,把残差压到局部峡谷底;返回 (x, r) 或 None。"""
    x = [float(x0[0]), float(x0[1])]
    r = rr(x[0], x[1])
    if r is None:
        return None
    v = _nrm(r)
    st = list(n_step0)
    st[0] = max(st[0]*max(x[0], 1.0), 1e-4)
    for _ in range(maxit):
        improved = False
        for ax in (0, 1):
            for sg in (1, -1):
                xp = list(x); xp[ax] += sg*st[ax]
                rp = rr(xp[0], xp[1])
                if rp is not None and _nrm(rp) < v*(1.0-1e-14):
                    x, r, v = xp, rp, _nrm(rp)
                    improved = True
                    break
            if improved:
                break
        if not improved:
            st[0] *= 0.5; st[1] *= 0.5
            if st[1] < 1e-8:
                break
    return (x, r)

def _refine_floor(rr, center, nstep_slide=(0.012, 0.008)):
    """对一个分支中心做: 模式下降(沿峡谷滑向底) + 牛顿;返回最优 (x,r) 或 None。"""
    best = None
    # 牛顿
    for seed in (center, (center[0]*1.01, center[1]+0.002),
                 (center[0]*0.99, center[1]-0.002)):
        hit = _newton(rr, seed)
        if hit is not None:
            return hit
    # 模式下降
    hit = _pattern_min(rr, center, n_step0=nstep_slide)
    if hit is not None and (best is None or _nrm(hit[1]) < _nrm(best[1])):
        best = hit
    return best

def solve_fixed_sc(s_c, p, nstep, prev=None,
                   n_coarse=(38, 26), n_basin=3, zoom_rounds=2):
    """在 (Ns0, psi_c) 平面解打靶问题。返回 dict(x, res, exact, norm) 或 None
    (全域找不到有限膜平衡)。策略:
      0) prev 延拓牛顿(扫描时沿解曲线免费精化, 命中精确根直接返回);
      1) 全域 log-Ns0 x lin-psi 粗采样, 按 |r| 排前 n_basin 个分支(去重);
      2) 每个分支做 zoom_rounds 轮局部放大网格(窗口指数收缩), 跟踪最低残差点;
      3) 每分支最低点做 牛顿->模式下降 抛光; 取全局最优;
      exact=True 当残差达到 ~1e-9(精确纯膜闭解); 否则 exact=False, res 为
      准闭解残差(峡谷底), 供与 COMSOL 对照(见文件头)。"""
    cache = {}
    psi_lo = math.pi/2.0 + 0.05
    psi_hi = math.pi/2.0 + 1.46
    rr = lambda a, b: (_safe_res(a, b, s_c, p, nstep, cache)
                       if psi_lo <= b <= psi_hi else None)
    Nsb = p*RC    # 特征膜力尺度(贴壁环向平衡 N_phi~-p*R_c, R_c=61.5 mm; 同量级替换旧参考圆估计)

    if prev is not None:
        for sd in (tuple(prev), (prev[0]*1.02, prev[1]+0.003),
                   (prev[0]*0.98, prev[1]-0.003)):
            hit = _newton(rr, sd)
            if hit is not None:
                return dict(x=hit[0], res=hit[1], exact=True,
                            norm=_nrm(hit[1]))

    # ---- 1) 全域粗采样 ----
    nNs, npsi = n_coarse
    Ns_lo, Ns_hi = 0.16*Nsb, 5.0*Nsb
    if prev is not None:
        Ns_lo = min(Ns_lo, 0.4*prev[0]); Ns_hi = max(Ns_hi, 2.0*prev[0])
    ns = [Ns_lo*(Ns_hi/Ns_lo)**(i/(nNs-1)) for i in range(nNs)]
    psid = [0.05 + 1.46*i/(npsi-1) for i in range(npsi)]   # psi = pi/2 + psid
    psis = [math.pi/2.0 + d for d in psid]
    pts = []
    for a in ns:
        for b in psis:
            r = rr(a, b)
            if r is not None:
                pts.append((_nrm(r), a, b, r))
    if not pts:
        return None
    pts.sort(key=lambda t: t[0])

    # ---- 分支去重(按 Ns0 相对距离聚类) ----
    basins = []
    for t in pts:
        nm, a, b, r = t
        for bk in basins:
            if abs(math.log(bk[0]/a)) < 0.35:
                break
        else:
            basins.append([a, b, nm])
        if len(basins) >= n_basin:
            break
    # ---- 2)+3) 逐分支放大 + 抛光 ----
    best = None
    for ba, bb, _nm in basins:
        cx = [ba, bb]
        nm = _nm
        r = None
        for _ in range(zoom_rounds):
            f = 0.30*max(1.0, cx[0]/Nsb)
            nsa = [cx[0]*(1.0-0.12), cx[0]*(1.0+0.12)]
            psa = [max(cx[1]-0.045, psi_lo), min(cx[1]+0.045, psi_hi)]
            if psa[1] - psa[0] < 2e-4:
                break
            n2, q2 = 18, 14
            loc = None
            for i in range(n2):
                aa = nsa[0]*(nsa[1]/nsa[0])**(i/(n2-1))
                for j in range(q2):
                    bbb = psa[0] + (psa[1]-psa[0])*j/(q2-1)
                    r2 = rr(aa, bbb)
                    if r2 is None:
                        continue
                    nm2 = _nrm(r2)
                    if loc is None or nm2 < loc[0]:
                        loc = (nm2, aa, bbb, r2)
            if loc is not None and loc[0] < nm:
                nm, cx[0], cx[1], r = loc[0], loc[1], loc[2], loc[3]
        hit = _refine_floor(rr, tuple(cx))
        if hit is not None and (best is None or _nrm(hit[1]) < _nrm(best[1])):
            best = hit
    if best is None:
        best = (tuple(cx), r)
    x, r = best
    return dict(x=list(x), res=r, exact=(_nrm(r) < 1e-9),
                norm=_nrm(r))

# ============================================================
#  4. 后处理与主流程
# ============================================================
def postprocess(out, p):
    S, r, z = out['S'], out['r'], out['z']
    Ns, lam1, lam2, q = out['Ns'], out['lam1'], out['lam2'], out['q']
    cs = out['sc_steps']
    curve = list(zip(r, z))
    poly = curve + [(ri, -zi) for ri, zi in reversed(curve)]
    V = 0.0
    for i in range(len(poly)-1):
        r1, z1 = poly[i]; r2, z2 = poly[i+1]
        V += (r1*r1 + r2*r2)*(z2 - z1)
    V = abs(0.5*math.pi*V)
    E = 0.0
    for i in range(len(S)-1):
        r0m = 0.5*(ref_r0_local(S[i]) + ref_r0_local(S[i+1]))
        l1m = 0.5*(lam1[i]+lam1[i+1]); l2m = 0.5*(lam2[i]+lam2[i+1])
        I1 = l1m**2 + l2m**2 + (1.0/(l1m*l2m))**2
        W = C10*(I1-3.0) + C20*(I1-3.0)**2 + C30*(I1-3.0)**3
        E += 2.0*math.pi*r0m*H0*W*(S[i+1]-S[i])
    E *= 2.0
    qq = q[1:cs+1]
    return dict(p=p, sc=out['sc_eff'], zc=max(z[:cs+1]) if cs > 0 else 0.0,
                Ns0=Ns[0], psi_c=out['psi_c'], z_end=out['res'][0],
                angle_res=out['angle_res'], cos_end=math.cos(out['psi_end']),
                q_peel=q[cs], qmin=min(qq) if qq else 0.0,
                psi_start=out['psi_start'], psi_end=out['psi_end'],
                r_start=out['r_start'], r_peel=out['r_peel'], r_end=out['r_end'],
                drds_free_start=out['drds_free_start'], rmax=max(r), E=E, V=V,
                lam1min=min(lam1), lam1max=max(lam1), Ns_min=min(Ns))

def fmt_row(d, exact):
    label = '精确闭解' if exact else '准闭解(残差见上)'
    return ('p=%7.1f Pa | s_c=%8.3f mm  z_c=%7.3f mm | Ns0=%8.2f N/m | psi_c=%6.3f rad | '
            'z_end=%+9.4f mm  angle_res=%+9.5f  cos_end=%+9.5f | q(sc)=%8.2f Pa | '
            'r_end=%7.2f mm | E_mem=%8.5f J  V=%7.3fe-5 m^3 | Ns_min=%8.2f N/m  %s'
            % (d['p'], d['sc']*1e3, d['zc']*1e3, d['Ns0'], d['psi_c'],
               d['z_end']*1e3, d['angle_res'], d['cos_end'], d['q_peel'], d['r_end']*1e3,
               d['E'], d['V']*1e5, d['Ns_min'], label))

def main():
    ap = argparse.ArgumentParser(description='WSL 充压段打靶(式2.5/3.3)')
    ap.add_argument('--p', type=float, default=400.0)
    ap.add_argument('--scan', type=str, default='')
    ap.add_argument('--scan-sc', type=str, default='')
    ap.add_argument('--fix-sc', type=float, default=None)
    ap.add_argument('--nstep', type=int, default=2400)
    ap.add_argument('--nstep-search', type=int, default=400)
    ap.add_argument('--out', type=str, default='')
    ap.add_argument('--quiet', action='store_true')
    a = ap.parse_args()
    mapping = local_material_mapping_check()

    def show_cfg():
        print('== 参数 ==')
        print('  Yeoh: C10=%g C20=%g C30=%g Pa | h0=%.2f mm | R_w=%.2f mm | R_c(贴壁中面)=%.2f mm'
              % (C10, C20, C30, H0*1e3, RW*1e3, RC*1e3))
        print('  中面参考环: 左 r=%.2f mm | 右 r=%.2f mm | 弧 R=%.2f mm | L0=%.3f mm'
              % (R_MID_L*1e3, R_MID_R*1e3, R_MID_A*1e3, L0*1e3))
        print('  COMSOL para=1 标定: p=%.1f Pa | s_c=%.3f mm (材料总长 %.3f mm; 当前接触长 %.3f mm)'
              % (P_COMSOL_REF, SC_COMSOL_REF*1e3, LC_MATERIAL_REF*1e3,
                 LC_CURRENT_COMSOL*1e3))
        print('  局部 S_hat=0 -> 全局 S=%.6f mm; 映射自检通过'
              % (S_CONTACT_CENTER_GLOBAL*1e3))
    show_cfg()

    rows = []
    def add_csv(d, sol):
        if a.out:
            rows.append((d, sol))
    def write_csv():
        if not a.out or not rows:
            return
        with open(a.out, 'w', newline='', encoding='utf-8') as fh:
            w = csv.writer(fh)
            w.writerow(['p_Pa','s_c_m','S_m','r0_m','z_m','r_m','psi_rad',
                        'Ns_Npm','lam1','lam2','q_Pa','seg'])
            for d, sol in rows:
                out = sol['out']
                for i, Sv in enumerate(out['S']):
                    seg = 'contact' if Sv <= out['sc_eff']+1e-9 else 'free'
                    w.writerow(['%.4f' % d['p'], '%.8f' % d['sc'],
                                '%.8f' % Sv, '%.8f' % ref_r0_local(Sv),
                                '%.8f' % out['z'][i], '%.8f' % out['r'][i],
                                '%.8f' % out['psi'][i], '%.8f' % out['Ns'][i],
                                '%.8f' % out['lam1'][i], '%.8f' % out['lam2'][i],
                                '%.8f' % out['q'][i], seg])
        print('构型 CSV 已写: %s' % a.out)

    def report(sol, scv, p, nstep_out):
        if sol is None:
            print('  s_c=%8.3f mm : 全域发散/无有限膜平衡(自由弧张力转压或不可反解)'
                  % (scv*1e3))
            return None, None
        x = sol['x']
        out = integrate_half(x[0], scv, x[1], p, nstep_out)
        d = postprocess(out, p)
        d['psi_c'] = x[1]
        exact = sol['exact'] and abs(d['z_end']) < 1e-4 and abs(d['angle_res']) < 1e-4
        if not exact:
            print('  s_c=%8.3f mm : 无精确纯膜闭解; 准闭解残差 |z_end|=%.2e m,'
                  ' |angle_res|=%.2e rad (需弯曲/CZM 补自由度, 见文件头)'
                  % (scv*1e3, abs(d['z_end']), abs(d['angle_res'])))
        print(' ' + fmt_row(d, exact))
        return d, dict(x=x, res=sol['res'], out=out)

    # ---- 沿 s_c 扫描(推荐默认用法) ----
    if a.scan_sc or (a.fix_sc is None and not a.scan):
        sc_list = ([float(x) for x in a.scan_sc.split(',') if x.strip()]
                   if a.scan_sc else [6e-3, 8e-3, 10.5e-3, 13e-3, 16e-3,
                                      20e-3, 25e-3, 30e-3])
        p = a.p
        print('== 沿贴壁弧半长 s_c 扫描(p=%g Pa; 搜索 nstep=%d, 出图 nstep=%d) =='
              % (p, a.nstep_search, a.nstep))
        print('   未知量 (Ns0, psi_c) 使 z(L0/2)=0 且 psi(L0/2)=3*pi/2; q(sc) 为诊断(见文件头)')
        prev = None
        nexact = 0
        for scv in sc_list:
            sol = solve_fixed_sc(scv, p, a.nstep_search, prev=prev)
            prev = None if sol is None else sol['x']
            d, payload = report(sol, scv, p, a.nstep)
            if d is not None and payload is not None:
                nexact += 1
                add_csv(d, payload)
        print('   找到 %d 个有限膜平衡(其中若干为精确闭解, 其余为准闭解;'
              ' 残差逐点打印)' % nexact)
        if not a.quiet:
            print('   - q(sc)>0: 剥离点仍需摩擦/CZM 锚定; q(sc)->0: 接近无黏附自然剥离点')
            print('   - 对照 COMSOL para=1: 取实测贴壁弧半长 s_c, 用 --fix-sc 单点精算')
        write_csv()
        return

    # ---- 单点 / 压力扫描(需 --fix-sc) ----
    ps = [float(x) for x in a.scan.split(',') if x.strip()] if a.scan else [a.p]
    if a.fix_sc is None:
        print('请给定贴壁弧半长: --fix-sc <m>(或使用默认 s_c 扫描模式)')
        return
    print('== 求解(固定 s_c=%.4f mm; 搜索 nstep=%d, 出图 nstep=%d) =='
          % (a.fix_sc*1e3, a.nstep_search, a.nstep))
    prev = None
    for p in ps:
        sol = solve_fixed_sc(a.fix_sc, p, a.nstep_search, prev=prev)
        prev = None if sol is None else sol['x']
        d, payload = report(sol, a.fix_sc, p, a.nstep)
        if d is not None and payload is not None:
            add_csv(d, payload)
    write_csv()

if __name__ == '__main__':
    main()
