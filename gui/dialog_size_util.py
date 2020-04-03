

class DialogSizeUtil:

    def __init__(self, q_app):
        self.app = q_app
        self.shape = self.app.desktop().screenGeometry()
        self.screen_width = self.shape.width()
        self.screen_height = self.shape.height()
        self.dialog_width = None
        self.dialog_height = None

    def get_window_width(self, width_percentage):
        self.dialog_width = round(self.screen_width * width_percentage)
        return self.dialog_width

    def get_window_height(self, height_percentage):
        self.dialog_height = round(self.screen_height * height_percentage)
        return self.dialog_height

    def center_dialog_width(self):
        return round(self.screen_width - self.dialog_width / 2)

    def center_dialog_height(self):
        return round(self.screen_height - self.dialog_height / 2)
