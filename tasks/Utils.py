from pathlib import Path
from tasks.History import HistoryHanlder
from tasks.Log import Log
from subprocess import run
import os

class Utils:
    def __init__(self, logSwitch = True) -> None:
        self.utilPath = Path(__file__)
        self.taskRoot = self.utilPath.parent.parent
        self.rRoot = self.taskRoot / 'R'
        self.meta = self.taskRoot / '.meta'
        self.logPath = self.taskRoot / 'log'
        self.historyPath = self.taskRoot / 'history'
        self.tmp = self.taskRoot / 'tmp'

        if logSwitch:
            self.log = Log(self.logPath, 'pipeline_log').logger
        else:
            self.log = Log(self.logPath, 'test_log').logger

        self.tmp.mkdir(parents=True, exist_ok=True)
        self.logPath.mkdir(parents=True, exist_ok=True)
        self.historyHanlder = HistoryHanlder(self.historyPath)

    def toAbsPath(self, pathStr) -> Path:
        path = Path(pathStr)
        if path.is_absolute(): 
            return path
        else:
            return self.taskRoot / pathStr
        
    def executeRtask(self, rtask):
        cmd_args = ['Rscript', rtask.script, *rtask.inputs, *rtask.outputs, *rtask.extraParam]
        cmd = [str(i) for i in cmd_args]
        if os.name == 'nt':
            cmd = [cmd[0]] + ["--max-mem-size=20000M"] + cmd[1:]
        executedCmd = ' '.join(cmd)
        self.log.info(f'executing cmd: {executedCmd}')
        self.process = run(cmd, capture_output=True)
        self.log.info(self.process.returncode)
        self.log.info(self.process.stdout)
        self.log.info(self.process.stderr)
        self.historyHanlder.putHistory(rtask, self.process)

        if self.process.returncode != 0:
            self.log.error(self.process.stderr)
            raise RuntimeError(self.process.stderr)


    