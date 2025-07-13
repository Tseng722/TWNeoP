import logging

class TWNeoP_Logger:
    def __init__(self, output_dir, name=__name__):
        self.output_dir = output_dir
        self.name = name
        self.logger = logging.getLogger(self.name)
        self.logger.setLevel(logging.DEBUG)

        if not self.logger.handlers:
            self._add_file_handler()

    def _add_file_handler(self):
        file_handler = logging.FileHandler(self.output_dir+"/log.log")
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)

        self.logger.addHandler(file_handler)

    def get_logger(self):
        return self
    
    def log_info(self,info_txt):
        self.logger.info(info_txt)
        print(info_txt)

    def log_error(self,info_txt):
        self.logger.error(info_txt)
        print(info_txt)