
var gulp = require('gulp');
var concat = require('gulp-concat');
var uglify = require('gulp-uglify');
var wrap = require('gulp-wrap');

var SOURCE_FILES = [
  'src/chromosome.js',
  'src/model-loader.js',
  'src/selector.js',
  'src/utils.js',
  'src/angular.adapter.js'
];

gulp.task('default', ['build-dev', 'release']);

gulp.task('release', function(){
  gulp.src(SOURCE_FILES)
    .pipe(concat('cyto-chromosome.min.js'))
    .pipe(wrap({src: 'build-wrap.template'}))
    .pipe(uglify())
    .pipe(gulp.dest('.'));
});

gulp.task('build-dev', function(){
  gulp.src(SOURCE_FILES)
    .pipe(concat('cyto-chromosome.js'))
    .pipe(wrap({src: 'build-wrap.template'}))
    .pipe(gulp.dest('.'));
});