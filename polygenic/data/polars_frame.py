"""A wrapper around a Polars DataFrame"""
import logging
from tqdm import tqdm
import sys
import os
import polars as pl
import pandas as pd
from enum import Enum

from logdecorator import log_on_start, log_on_end, log_on_error

class GwasSchema(Enum):

      """GwasSchema v0.1.0"""

      def __new__(cls, *args, **kwds):
          value = len(cls.__members__) + 1
          obj = object.__new__(cls)
          obj._value_ = value
          return obj
      def __init__(self, colname, aliases, cast):
          self.colname = colname
          self.aliases = aliases
          self.cast = cast

      CHROMOSOME = "chromosome", ["chromosome", "chrom", "chr"], pl.Utf8
      POSITION = "position", ["position", "pos"], pl.Int64
      GNOMAD_ID = "gnomad_id", ["gnomad_id"], pl.Utf8
      RSID = "rsid", ["rsid", "rs"], pl.Utf8
      REFERENCE = "reference", ["reference", "ref"], pl.Utf8
      ALTERNATIVE = "alternative", ["alternative", "alt"], pl.Utf8
      EFFECT_ALLELE = "effect_allele", ["effect_allele", "effect", "alt"], pl.Utf8
      NON_EFFECT_ALLELE = "non_effect_allele", ["non_effect_allele"], pl.Utf8
      PVALUE = "pvalue", ["pvalue", "pval", "p"], pl.Float64
      PVALUE_HETEROGENEITY = "pvalue_heterogeneity", ["pvalue_heterogeneity", "pval_heterogeneity", "p_heterogeneity"], pl.Float64
      BETA = "beta", ["beta"], pl.Float64
      SE = "se", ["se"], pl.Float64
      ZSCORE = "zscore", ["zscore", "z"], pl.Float64
      OR = "or", ["or"], pl.Float64
      INFO = "info", ["info"], pl.Float64
      MAF = "maf", ["maf"], pl.Float64
      AF = "af", ["af"], pl.Float64
      LOW_CONFIDENCE = "low_confidence", ["low_confidence", "low_conf"], pl.Boolean
    
class PolarsFrame:

  @log_on_start(logging.DEBUG, "Initializing PolarsFrame from {csv}")
  def __init__(self, csv, delimiter='\t'):
    self.df = self.readcsv(csv, delimiter)

  @log_on_start(logging.DEBUG, "Getting dataframe")
  def get_dataframe(self):
    return self.df

  @log_on_start(logging.DEBUG, "Reading {csv}")
  @log_on_end(logging.DEBUG, "Read {csv}")
  def readcsv(self, csv, delimiter='\t'):
    CHUNK_SIZE = 10**5
    CHUNK_ROWS = 500
    initial_chunk = pd.read_csv(filepath_or_buffer = csv, sep = delimiter, nrows = CHUNK_ROWS, dtype=str)
    initial_chunk_size = len(initial_chunk.to_csv(index=False))
    output_df = pl.from_pandas(initial_chunk)
    total_rows = CHUNK_ROWS * int(os.path.getsize(csv) / initial_chunk_size * CHUNK_ROWS * 2.5 / CHUNK_SIZE) + 1
    with tqdm(total = total_rows, file = sys.stdout, leave=False) as pbar:
      for i,chunk in enumerate(pd.read_csv(csv, sep = delimiter, chunksize=CHUNK_SIZE, low_memory=False, dtype=str)):
        output_df.extend(pl.from_pandas(chunk))
        pbar.set_description('Reading csv chunks (estimated): %d' % ((1 + i) * CHUNK_ROWS))
        pbar.update(CHUNK_ROWS)
    return output_df

  @log_on_start(logging.DEBUG, "Getting column names")
  def get_cols(self):
    return self.df.columns

  def get_col_by_aliases_list(self, aliases, keyword = None):
    for alias in aliases:
      col_name = self.get_col_by_alias(alias, keyword)
      if col_name is not None:
        return col_name
    return 'empty'

  def get_col_by_alias(self, alias, keyword = None):
    if alias is None:
      return None
    alias_variants = [alias, alias.lower(), alias.upper(), alias.capitalize()]
    if keyword is not None:
      column_names = [column for column in self.df.columns if keyword in column] + [column for column in self.df.columns if keyword not in column]
    else:
      column_names = self.df.columns
    for alias_variant in alias_variants:
      if alias_variant in column_names:
        return alias_variant
    for alias_variant in alias_variants:
      for col in self.df.columns:
        if col.startswith(alias_variant):
          return alias_variant
    return None

  @log_on_start(logging.DEBUG, "Converting to gwas object")
  def convert_to_gwas(self, col_mappings: dict):
    self.df = self.df.with_column(pl.lit(None).alias('empty')).lazy().select([
        pl.col(self.get_col_by_aliases_list([col_mappings.get(col.colname)] + col.aliases)).cast(col.cast).alias(col.colname) for col in GwasSchema
      ]).collect()

