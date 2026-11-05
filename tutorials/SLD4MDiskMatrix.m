classdef SLD4MDiskMatrix < handle
% Lossless double-precision storage adapter for the unchanged SLD4M.m.
% The binary file contains rows of the original matrix, in original order.
% Only the array operations used by the author's default example are exposed.
  properties (SetAccess = private)
    filename; nrows; filecols; columns; column_sums;
  end
  properties (Access = private)
    fid; cached_rows = []; cached_data = [];
  end
  methods
    function obj = SLD4MDiskMatrix(filename, nrows, ncols, column_sums, columns)
      if nargin < 5, columns = 1:ncols; end
      obj.filename = filename; obj.nrows = double(nrows); obj.filecols = double(ncols);
      obj.columns = double(columns); obj.column_sums = double(column_sums);
      obj.fid = fopen(filename, 'rb', 'ieee-le');
      if obj.fid < 0, error('Cannot open reference matrix: %s', filename); end
    end
    function delete(obj)
      if ~isempty(obj.fid) && obj.fid >= 0, fclose(obj.fid); end
    end
    function varargout = size(obj, dim)
      s = [obj.nrows, numel(obj.columns)];
      if nargin > 1
        if dim <= 2, varargout{1} = s(dim); else, varargout{1} = 1; end
      elseif nargout <= 1, varargout{1} = s;
      else, varargout{1} = s(1); varargout{2} = s(2);
      end
    end
    function n = end(obj, k, nd)
      if nd ~= 2, error('Only two-dimensional indexing is supported'); end
      n = size(obj, k);
    end
    function y = subsref(obj, s)
      if ~strcmp(s(1).type, '()')
        y = builtin('subsref', obj, s); return;
      end
      if numel(s(1).subs) ~= 2, error('Use two-dimensional indexing'); end
      rows = s(1).subs{1}; cols = s(1).subs{2};
      if ischar(cols) && strcmp(cols, ':'), cols = 1:numel(obj.columns); end
      selected = obj.columns(cols);
      if ischar(rows) && strcmp(rows, ':')
        % Native code uses this view only with mean(), for enrichment totals.
        y = SLD4MDiskMatrix(obj.filename, obj.nrows, obj.filecols, obj.column_sums, selected);
      else
        rows = rows(:);
        if any(rows < 1 | rows > obj.nrows | rows ~= fix(rows))
          error('Reference row index out of range');
        end
        if ~isequal(rows, obj.cached_rows)
          lo = min(rows); hi = max(rows);
          if fseek(obj.fid, (lo-1)*obj.filecols*8, 'bof') ~= 0, error('Seek failed'); end
          [block, count] = fread(obj.fid, [obj.filecols, hi-lo+1], '*double');
          if count ~= obj.filecols*(hi-lo+1), error('Truncated reference matrix'); end
          obj.cached_data = block(:, rows-lo+1)'; obj.cached_rows = rows;
        end
        y = ones(numel(rows), numel(selected));
        nonconstant = selected > 0;
        y(:,nonconstant) = obj.cached_data(:,selected(nonconstant));
      end
      if numel(s) > 1, y = subsref(y, s(2:end)); end
    end
    function y = sum(obj, dim)
      if nargin > 1 && dim ~= 1, error('Only column sums are supported'); end
      y = repmat(obj.nrows, 1, numel(obj.columns));
      realcols = obj.columns > 0;
      y(realcols) = obj.column_sums(obj.columns(realcols));
    end
    function y = mean(obj, dim)
      if nargin > 1 && dim ~= 1, error('Only column means are supported'); end
      y = sum(obj) / obj.nrows;
    end
    function y = horzcat(prefix, obj)
      if ~isa(obj, 'SLD4MDiskMatrix') || ~isnumeric(prefix) || ...
          ~isequal(size(prefix), [obj.nrows, 1]) || any(prefix ~= 1)
        error('Only the native all-ones intercept prefix is supported');
      end
      y = SLD4MDiskMatrix(obj.filename, obj.nrows, obj.filecols, ...
          obj.column_sums, [0, obj.columns]);
    end
  end
end
