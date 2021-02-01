# -*- coding: utf-8 -*-
#!/usr/bin/env python
# Author: Ziwei Xue
#
# ------------------------- #
# Python Modules
# ------------------------- #

import sqlite3

def sqlite_connect(db):
	return sqlite3.connect(db)

def sqlite_delete_gtf_table(conn,tbname):
	c = conn.cursor()
	c.execute("""DROP TABLE {} """.format(tbname))
	conn.commit()
	conn.close()

def sqlite_create_gtf_table(conn,tbname):
	gtf_format = """ (seqname VARCHAR(10) NOT NULL,
                source VARCHAR(20) NOT NULL,
                feature VARCHAR(20) NOT NULL,
                start INTEGER NOT NULL,
                end INTEGER NOT NULL,
                score VARCHAR(10) NOT NULL,
                strand VARCHAR(1) NOT NULL,
                frame VARCHAR(1) NOT NULL,
                gene_id VARCHAR(25),
                gene_version VARCHAR(10),
                transcript_id VARCHAR(20),
                transcript_version VARCHAR(20),
                transcript_name VARCHAR(20),
                transcript_source VARCHAR(20),
                transcript_biotype VARCHAR(20),
                transcript_support_level VARCHAR(20),
                exon_number VARCHAR(2),
                exon_id VARCHAR(20),
                exon_version VARCHAR(20),
                exon_name VARCHAR(20),
                exon_source VARCHAR(20),
                exon_biotype VARCHAR(20),
                protein_id VARCHAR(20),
                protein_version VARCHAR(20),
                protein_name VARCHAR(20),
                protein_source VARCHAR(20),
                protein_biotype VARCHAR(20),
                ccds_id VARCHAR(20),
                ccds_version VARCHAR(20),
                ccds_name VARCHAR(20),
                ccds_source VARCHAR(20),
                ccds_biotype VARCHAR(20),
                gene_name VARCHAR(25),
                gene_source VARCHAR(10),
                gene_biotype VARCHAR(50),
                tag VARCHAR(10))"""
	c = conn.cursor()
	c.execute("""CREATE TABLE {} """.format(tbname) + gtf_format)
	conn.commit()
	conn.close()


def sqlite_insert_gtf_table(conn, tabname, data):
	c = conn.cursor()
	for i in data:
		const, attr = i[0], i[1]
		sql = "INSERT INTO {}(seqname, source, feature, start, end, score, strand, frame, ".format(tabname) + ','.join(attr.keys()) + ") VALUES('%s', '%s', '%s','%d', %d, '%s', '%s', '%s', " + ",".join(["%s"] * len(attr)) + ")"
		c.execute(sql % tuple(const[:3] + [int(const[3]) ,int(const[4])] + const[5:]  + list(map(lambda x:"'{}'".format(x), attr.values()))))
	conn.commit()
	conn.close()

def sqlite_query_gtf_by_genename(conn, tabname, genename):
	c = conn.cursor()
	c.execute("SELECT seqname, source, feature, start, end, score, strand, frame, gene_name FROM {} WHERE gene_name='{}'".format(tabname, genename))
	result = c.fetchall()
	return list(map(lambda x:list(x)[:8] + [{"gene_name":x[8]}], result))               

def sqlite_query_gtf_by_region(conn, tabname, chrom, start, end):
	c = conn.cursor()
	c.execute("SELECT seqname, source, feature, start, end, score, strand, frame, gene_name FROM {} WHERE seqname='{}' AND start>{} AND end<{}".format(tabname, chrom, start, end))
	result = c.fetchall()
	return list(map(lambda x:list(x)[:8] + [{"gene_name":x[8]}], result))     

def sqlite_query_gtf_by_ensembl_geneid(conn, tabname, geneid):
	c = conn.cursor()
	c.execute("SELECT seqname, source, feature, start, end, score, strand, frame, gene_name FROM {} WHERE gene_id='{}'".format(tabname, genename))
	result = c.fetchall()
	return list(map(lambda x:list(x)[:8] + [{"gene_name":x[8]}], result))    