import time


class SequenceAnimator(object):
    def __init__(self, time, map_meshes, color_maps):
        self.angle = 0
        self.animation_state = False
        self.start_time = 0
        self.time = time
        self.map_meshes = map_meshes
        self.time_flag = False
        self.pos = 0
        self.color_maps = color_maps

    def __call__(self, bg_plotter, **kwargs):
        if self.animation_state:
            if int(self.start_time+time.time()) % self.time == 0:
                if self.time_flag:
                    bg_plotter.clear()
                    self.pos += 1
                    map_pos = self.pos%len(self.map_meshes)
                    bg_plotter.add_mesh(self.map_meshes[map_pos], cmap=self.color_maps[map_pos], **kwargs)
                    self.time_flag = False
            else:
                self.time_flag = True
        else:
            self.start_time = time.time()