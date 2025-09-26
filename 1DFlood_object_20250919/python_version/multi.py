# 各数値計算方は，fluxを出すためだけに変更する
class Euler:
    
    # 水深の時間発展
    def update_depth(self,Element,dt):
        uppdated_depth = Element.get_variable_depth()+Element.calc_increment()*dt
        Element.set_depth(uppdated_depth)