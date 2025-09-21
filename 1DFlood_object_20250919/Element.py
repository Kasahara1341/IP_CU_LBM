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
    def set_AdamsBashforth4th(self):
        self.time_evo = AdamsBashforth4th()

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
    def __init__(self) -> None:
        # self.prev_time = [0,0,0,0]    # 可変dtに対応するなら必要
        self.coeff = [55.0/24.0, -59.0/24.0, 37.0/24.0, -9.0/24.0]
        self.prev_dh = [0,0,0,0]
    # incrementリストを更新 [0]が最新，[3]が一番古い
    def update_increment_list(self,Element):
        for i in range(3):
            self.prev_dh[i+1] = self.prev_dh[i]
        dh = Element.calc_increment()
        self.prev_dh[0] = dh
    # 水深の時間発展
    def update_depth(self,Element,dt):
        self.update_increment_list(Element)
        sum_increment = 0
        for i in range(4):
            sum_increment += self.coeff[i]*self.prev_dh[i]
        uppdated_depth = Element.get_variable_depth() + dt*sum_increment
        Element.set_depth(uppdated_depth)