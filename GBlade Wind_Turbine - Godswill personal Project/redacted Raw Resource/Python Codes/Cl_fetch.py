import sqlite3
conn=sqlite3.connect("airfoil.db")
cur=conn.cursor()
#cur.execute("create table if not exists Apha_vs_Cl (Apha integer primary key, Cl)")
conn.commit()
def Cl_search(Al):
    try:
        cur.execute("SELECT Apha FROM Apha_vs_Cl")
        xlist = cur.fetchall()
        x=[]
        for tup in xlist:
            x1=float(str(tup).replace('(','').replace(',)',''))
            x.append(x1)
        Cl=x[min(range (len(x)), key = lambda i: abs(x[i]-Al))]
        cur.execute("SELECT* FROM Apha_vs_Cl WHERE Apha=?", (Cl, ))
        Clist = cur.fetchall()
        Cl=Clist[0]
        return Cl[0]
    except:
        conn.rollback()
def close_db():
    conn.close()
