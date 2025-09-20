import numpy as np
class Node:
    def __init__(self):
        self.q = None
        self.down_element =None
        self.up_element = None
        

    def set_upelement(self,up_element):
        self.up_element = up_element

    def set_downelement(self,down_element):
        self.down_element = down_element

    def get_upelement(self):
        return self.up_element
    
    def get_downelement(self):
        return self.down_element

    def set_q(self,value):
        self.q = value

    def get_variable_q(self):
        return self.q
    
    # Nodeにおけるフラックス
    def calc_flux(self):
        #dx = abs(self.down_element.get_position()-self.up_element.get_position())
        flux = self.q/150
        return flux
    
    # 運動方程式(拡散波近似)
    def solve_momentum_equation(self):
        Hup     = self.up_element.get_elev() + self.up_element.get_variable_depth()
        Hdn     = self.down_element.get_elev() + self.down_element.get_variable_depth()
        hup     = self.up_element.get_variable_depth()
        hdn     = self.down_element.get_variable_depth()
        Bup = self.up_element.get_width()
        Bdn = self.down_element.get_width()
        
        # h = (hup+hdn)/2.0 original
        # B = (Bup+Bdn)/2.0 original
        h = (hup) # 0919 23:52
        B = (Bup) # 0919 23:52

        h = np.max([0.0,h])    # 水位が負値をとったときの処理
        n = (self.up_element.get_coeff()+self.down_element.get_coeff())/2.0
        dx    = abs(self.down_element.get_position()-self.up_element.get_position()) 
        grad = np.abs((Hup-Hdn)/dx)
        self.q = B*(1/n)*h**(5/3)*np.sqrt(grad)*np.sign((Hup-Hdn)/dx)
