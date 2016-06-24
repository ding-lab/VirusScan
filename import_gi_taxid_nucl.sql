\! echo 'Loading gi_taxid_nucl.dmp - this can take a very long time'
DROP TABLE IF EXISTS `gi_taxid_nucl`;
CREATE TABLE `gi_taxid_nucl` (
  `gi_taxid_nucl_id` int(10) unsigned NOT NULL auto_increment,
  `gi` int(10) unsigned default NULL,
  `tax_id` int(10) unsigned default NULL,
  PRIMARY KEY  (`gi_taxid_nucl_id`),
  KEY `tax_id` (`tax_id`),
  KEY `gi` (`gi`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


LOAD DATA LOCAL INFILE '/directory/gi_taxid_nucl.dmp'
    INTO TABLE gi_taxid_nucl
    FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n'
    (gi,tax_id);
