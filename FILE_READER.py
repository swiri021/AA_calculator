__author__ = 'junheeyun'
#import MySQLdb
#from sys import argv
#import subprocess
#import os

class FILE_READ:

	def returning(self):
		return self.result_arr

	def returning_selected(self,sel1,sel2):
		temp_arr = []
		for a in self.result_arr:
			temp_arr.append([a[sel1], a[sel2]])
		return temp_arr



	def __init__(self, file_name, spliter):
		self.result_arr = []
		f = open(file_name, "r")

		for x in f.readlines():
			temp = x.split(spliter)
			cleared =[]
			for y in temp:
				y = y.strip()
				cleared.append(y)

			self.result_arr.append(cleared)

		f.close()

	def __del__(self):
		print "terminated file_read"
