from RK6th import Runge_Kutta_6th
import numpy as np
class Element: #格子のオブジェクト
    def __init__(self, position,elev,coeff,width): # コンストラクタ(selfがオブジェクトの個性的な感じ，それ以降の変数は個性の内訳)
        self.position = position # 格子点の位置
        self.elev = elev # 格子点の標高
        self.coeff = coeff  # 粗度係数
        self.width = width #川幅
        self.depth = None  # Noneのやつはあとでdefで決める必要あり
        self.old_depth = None
        self.upnoads = []
        self.dnnoads = []
        self.time_evo  = None

    #数値計算法の選択
    def set_Euler(self):
        self.time_evo = Euler()
    #数値計算法の選択
    def set_AdamsBashforth4th(self,dt):
        self.time_evo = AdamsBashforth4th(dt)
    def set_Runge_Kutta_4th(self):
        self.time_evo = Runge_Kutta_4th()        
    def set_Runge_Kutta_6th(self):
        self.time_evo = Runge_Kutta_6th()        

    # 質量保存則
    def solve_mass_equation(self,dt):
        self.time_evo.update_depth(self,dt)   # hの更新
    
    def calc_increment(self):
        fluxin = 0.0
        fluxout = 0.0

        for upnoad in self.upnoads:
            fluxin += upnoad.calc_flux()
        
        for dnnoad in self.dnnoads:
            fluxout += dnnoad.calc_flux()
            
        dh = (fluxin - fluxout)/self.width
        return dh

    #def set_upnoad(self,upnoad):
    #    self.upnoad = upnoad

    def set_upnoad(self,upnoad):
        self.upnoads.append(upnoad)
        
    def set_dnnoad(self,dnnoad):
        self.dnnoads.append(dnnoad)
    
    def set_depth(self,value):
        self.depth = value
    
    def set_old_depth(self,value):   # 現時刻の水深
        self.old_depth = value

    def get_variable_depth(self):
        return self.depth

    def get_variable_old_depth(self):
        return self.old_depth
    
    def get_position(self):
        return self.position
    
    def get_elev(self):
        return self.elev
    
    def get_width(self):
        return self.width
    
    def get_coeff(self):
        return self.coeff


class Euler:
    
    # 水深の時間発展
    def update_depth(self,Element,dt):
        uppdated_depth = Element.get_variable_depth()+Element.calc_increment()*dt
        Element.set_depth(uppdated_depth)

class AdamsBashforth4th:
    def __init__(self,dt) :
        self.prev_time = [-1*dt,-2*dt,-3*dt,-4*dt]    # 可変dtに対応するなら必要
        self.coeff = [55.0/24.0, -59.0/24.0, 37.0/24.0, -9.0/24.0]
        self.prev_dh = [0,0,0,0]
    # Lagrange補間多項式の積分をによる係数計算
    def integrated_Lagrange4(self,index,dt):
        t0 = self.prev_time[(0+index)%4]
        t1 = self.prev_time[(1+index)%4]
        t2 = self.prev_time[(2+index)%4]
        t3 = self.prev_time[(3+index)%4]
        time = max(t0,t1,t2,t3) + dt
        coeff4 = 3.0*time**4.0 - 4.0*(t1+t2+t3)*time**3.0 + \
            6.0*(t1*t2 + t2*t3 + t3*t1)*time**2.0 - 12.0*t1*t2*t3*time
        coeff4 /= 12.0*(t0 - t1)*(t0 - t2)*(t0 - t3)
        return coeff4
    # 現在のprev_timeに基づいて係数を更新
    def reset_coeff(self,dt):
        sum_coeff = 0 
        for i in range(4):
            self.coeff[i] = self.integrated_Lagrange4(i,dt) - self.integrated_Lagrange4(i,0)
            self.coeff[i] /= dt
        sum_coeff = np.sum(self.coeff)
        # 係数の和が1となるように正規化
        for i in range(4):
            self.coeff[i] /= sum_coeff
    # incrementリストを更新 [0]が最新，[3]が一番古い
    def update_increment_list(self,Element,dt):
        for i in range(3):
            self.prev_dh[3-i] = self.prev_dh[2-i]
            self.prev_time[3-i] = self.prev_time[2-i]
        dh = Element.calc_increment()
        self.prev_dh[0] = dh
        self.prev_time[0] = self.prev_time[0] + dt
    # 水深の時間発展
    def update_depth(self,Element,dt):
        self.update_increment_list(Element,dt)
        # self.reset_coeff(dt)
        sum_increment = 0
        for i in range(4):
            sum_increment += self.coeff[i]*self.prev_dh[i]
        uppdated_depth = Element.get_variable_depth() + dt*sum_increment
        Element.set_depth(uppdated_depth)

class Runge_Kutta_4th(object):
    def __init__(self):
        self.stage = 0
        self.hold = None
        self.increments = []

    def update_stage(self):
        if self.stage < 3:
            self.stage += 1
        else:
            self.stage = 0
    
    def set_hold(self,value):
        self.hold = value
    
    def update_stage0_variables(self,Element,dt):
        k1 = Element.calc_increment()
        uppdated_depth = self.hold + 0.5*k1*dt
        Element.set_depth(uppdated_depth)
        self.increments.append(k1)
        self.update_stage()
    
    def update_stage1_variables(self,Element,dt):
        k2 = Element.calc_increment()
        uppdated_depth = self.hold + 0.5*k2*dt
        Element.set_depth(uppdated_depth)
        self.increments.append(k2)
        self.update_stage()
    
    def update_stage2_variables(self,Element,dt):
        k3 = Element.calc_increment()
        uppdated_depth = self.hold + k3*dt
        Element.set_depth(uppdated_depth)
        self.increments.append(k3)
        self.update_stage()

    def update_stage3_variables(self,Element,dt):
        k4 = Element.calc_increment()
        self.increments.append(k4)
        self.update_stage()

    def update_depth(self,Element,dt):
        if self.stage == 0:
            self.hold = Element.get_variable_depth()
            self.update_stage0_variables(Element,dt)

        elif self.stage==1:
            self.update_stage1_variables(Element,dt)
        
        elif self.stage==2:
            self.update_stage2_variables(Element,dt)

        else:
            self.update_stage3_variables(Element,dt)
            uppdated_depth = self.hold + (1/6)*(self.increments[0]+2*self.increments[1]+2*self.increments[2]+self.increments[3])*dt
            Element.set_depth(uppdated_depth)
            self.increments = []
        