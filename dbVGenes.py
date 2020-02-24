__author__ = 'wilsonp'
import sqlite3 as db
# first need connect to a database
conn = db.connect('test.db')

#  then need to create a cursor that lets you traverse the database

cursor = conn.cursor()
conn = db.connect('test.db')

#  then need to create a cursor that lets you traverse the database


cursor.execute('create table films(title text, year text, director text)')

cursor.execute('drop table if exists temps')
cursor.execute('create table temps(dat text, temp int)')
cursor.execute('insert into temps values("12/1/2011", 35)')
cursor.execute('insert into temps values("12/2/2011", 13)')
cursor.execute('insert into temps values("12/3/2011", 40)')
cursor.execute('insert into temps values("12/4/2011", 28)')
cursor.execute('insert into temps values("12/5/2011", 30)')
cursor.execute('insert into temps values("12/6/2011", 34)')
conn.commit()

conn.row_factory = db.Row
cursor.execute('select * from temps')
rows = cursor.fetchall()
for row in rows:
    print('%s %s' % (row[0], row[1]))

#can use sql to determine averages
cursor.execute('select avg(temp) from temps')
row = cursor.fetchone()
print('The average temp for the week was: %s' % row[0])
cursor.execute('delete from temps where temp = 40')

# need to requery and fetch data to display with changes
cursor.execute('select * from temps')
rows = cursor.fetchall()
for row in rows:
    print('%s %s' % (row[0], row[1]))

