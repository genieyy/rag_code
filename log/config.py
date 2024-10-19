import logging  
import json  

class JsonFormatter(logging.Formatter):  
    def format(self, record):  
        # 构建日志记录字典  
        log_record = {  
            "time": self.formatTime(record, self.datefmt), 
            "level": record.levelname,
            "name": record.name,  
            "message": record.msg,  
        }  
        
        if record.exc_info:  
            log_record["exc_info"] = self.formatException(record.exc_info)  
        
        return json.dumps(log_record)  

def setup_logging(record_path:str):  
    # 创建自定义JSON格式化器  
    json_formatter = JsonFormatter()  

    # 创建日志处理器  
    file_handler = logging.FileHandler(record_path)  
    file_handler.setFormatter(json_formatter)  

    console_handler = logging.StreamHandler()  
    console_handler.setFormatter(json_formatter)  

    # 获取根logger并添加处理器  
    root_logger = logging.getLogger()  
    root_logger.setLevel(logging.INFO)  
    root_logger.addHandler(file_handler)  
    root_logger.addHandler(console_handler)  
    return root_logger


