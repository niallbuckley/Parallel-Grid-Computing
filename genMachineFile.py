#!/usr/bin/env python
import threading, time, socket

prefix = "csg21-"
start = 2
end = 48
filename = "machines"

mFile = open(filename, "w")
lock = threading.Lock()
count = 0
threads = []

class ThreadChecker(threading.Thread):
	def __init__(self, hostNum):
		threading.Thread.__init__(self)
		self.hostNum = hostNum

	def run(self):
		global count
		hostName = prefix + str(self.hostNum).zfill(2) + ".ucc.ie"

		s = socket.socket()
		try:
			s.connect((hostName, 22))
			s.close()
			with lock:
				mFile.writelines(hostName+"\n")
				count += 1
		except:
			pass

if __name__ == '__main__':
	for i in range(start, end):
		thread = ThreadChecker(i)
		threads.append(thread)
		thread.start()

	for t in threads:
		t.join()

	mFile.close()
	print "%d machine(s) out of %d found"%(count,end - start)